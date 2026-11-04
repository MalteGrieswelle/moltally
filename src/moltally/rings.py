"""Enumerate and count the ring systems (fused rings merged) of a SMILES file."""

import argparse
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger

from moltally.common import (
    CHUNK,
    PROGRESS_EVERY,
    DATA,
    draw_top_structures,
    parallel_map,
    plot_distribution,
    read_chunks,
    write_counts,
    write_top_smiles,
)


def ring_systems(mol: Chem.Mol) -> list[tuple[int, ...]]:
    """Atom indices of each ring system.

    Rings sharing two or more atoms (fused or bridged) form one system; rings
    that share a single atom (spiro) or none stay separate systems.
    """
    systems: list[set[int]] = []
    for ring in mol.GetRingInfo().AtomRings():
        merged = set(ring)
        rest = []
        for system in systems:
            if len(system & merged) >= 2:
                merged |= system
            else:
                rest.append(system)
        systems = rest + [merged]
    return [tuple(sorted(system)) for system in systems]


def ring_smiles(mol: Chem.Mol, atom_ids: tuple[int, ...]) -> str | None:
    """Canonical SMILES of one ring system, stripped of substituents.

    Keeps the ring atoms plus atoms attached through exocyclic double bonds
    (C=O, C=S, ...). Hydrogens are re-assigned so that removed substituents
    become hydrogens, e.g. N-methylpyrrole -> c1cc[nH]c1. Returns None if the
    fragment cannot be sanitized.
    """
    ring = set(atom_ids)
    keep = set(ring)
    for i in ring:
        for bond in mol.GetAtomWithIdx(i).GetBonds():
            other = bond.GetOtherAtomIdx(i)
            if other not in ring and bond.GetBondType() == Chem.BondType.DOUBLE:
                keep.add(other)

    rw = Chem.RWMol(mol)
    for idx in keep:
        atom = rw.GetAtomWithIdx(idx)
        removed = sum(1 for n in atom.GetNeighbors() if n.GetIdx() not in keep)
        atom.SetNumExplicitHs(atom.GetTotalNumHs() + removed)
        atom.SetNoImplicit(True)
        atom.SetIsotope(0)
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    for idx in sorted(set(range(rw.GetNumAtoms())) - keep, reverse=True):
        rw.RemoveAtom(idx)
    for bond in rw.GetBonds():
        bond.SetStereo(Chem.BondStereo.STEREONONE)
    try:
        frag = rw.GetMol()
        Chem.SanitizeMol(frag)
    except Exception:
        return None
    return Chem.MolToSmiles(frag, isomericSmiles=False)


def _count_chunk(lines: list[str]) -> tuple[Counter, Counter, int]:
    RDLogger.DisableLog("rdApp.*")
    occurrences, molecules, failed = Counter(), Counter(), 0
    for line in lines:
        mol = Chem.MolFromSmiles(line.split()[0])
        if mol is None:
            continue
        smis = []
        for system in ring_systems(mol):
            smi = ring_smiles(mol, system)
            if smi is None:
                failed += 1
            else:
                smis.append(smi)
        occurrences.update(smis)
        molecules.update(set(smis))
    return occurrences, molecules, failed


def count_rings(path: Path, workers: int | None = None):
    occurrences, molecules, failed = Counter(), Counter(), 0
    for i, (occ, mol, fail) in enumerate(parallel_map(_count_chunk, read_chunks(path), workers), 1):
        occurrences.update(occ)
        molecules.update(mol)
        failed += fail
        if i % PROGRESS_EVERY == 0:
            print(f"{i * CHUNK:>11,} molecules, {len(occurrences)} ring systems", flush=True)
    return occurrences, molecules, failed


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, default=DATA / "chembl37.smi")
    p.add_argument("--outdir", type=Path, default=Path("results"))
    p.add_argument("--top", type=int, default=100)
    p.add_argument("--workers", type=int, default=None)
    a = p.parse_args()

    occ, mols, failed = count_rings(a.input, a.workers)
    print(f"{len(occ)} distinct ring systems; {failed} ring system instances could not be sanitized")
    a.outdir.mkdir(parents=True, exist_ok=True)
    stem = a.input.stem
    write_counts(a.outdir / f"{stem}_rings.tsv", occ, mols)
    plot_distribution(occ, mols, a.outdir / f"{stem}_ring_frequency.png", a.top, "ring systems")
    draw_top_structures(occ, mols, a.outdir / f"{stem}_top{a.top}_rings.png", a.top)
    write_top_smiles(a.outdir / "ring.smi", occ, mols, a.top)


if __name__ == "__main__":
    main()
