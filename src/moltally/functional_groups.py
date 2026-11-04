"""Count functional groups (Ertl IFG) over a SMILES file."""

import argparse
from collections import Counter
from functools import lru_cache, partial
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Contrib.IFG.ifg import identify_functional_groups

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


@lru_cache(maxsize=None)
def normalize_fg(smiles: str) -> str:
    """Canonical, aromaticity-, isotope- and stereo-free form of an IFG type string."""
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    for atom in mol.GetAtoms():
        atom.SetIsAromatic(False)
        atom.SetIsotope(0)
    for bond in mol.GetBonds():
        bond.SetIsAromatic(False)
        if bond.GetBondType() == Chem.BondType.AROMATIC:
            bond.SetBondType(Chem.BondType.SINGLE)
    mol.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(mol)
    return Chem.MolToSmiles(mol, isomericSmiles=False)


def _count_chunk(lines: list[str], exclude_ring: bool) -> tuple[Counter, Counter]:
    RDLogger.DisableLog("rdApp.*")
    occurrences, molecules = Counter(), Counter()
    for line in lines:
        mol = Chem.MolFromSmiles(line.split()[0])
        if mol is None:
            continue
        groups = identify_functional_groups(mol)
        if exclude_ring:
            # drop functional groups with any atom that is part of a ring
            ring_info = mol.GetRingInfo()
            groups = [
                fg
                for fg in groups
                if not any(ring_info.NumAtomRings(i) for i in fg.atomIds)
            ]
        types = [normalize_fg(fg.type) for fg in groups]
        occurrences.update(types)
        molecules.update(set(types))
    return occurrences, molecules


def count_functional_groups(
    path: Path, workers: int | None = None, exclude_ring: bool = True
) -> tuple[Counter, Counter]:
    occurrences, molecules = Counter(), Counter()
    worker = partial(_count_chunk, exclude_ring=exclude_ring)
    for i, (occ, mol) in enumerate(parallel_map(worker, read_chunks(path), workers), 1):
        occurrences.update(occ)
        molecules.update(mol)
        if i % PROGRESS_EVERY == 0:
            print(f"{i * CHUNK:>11,} molecules, {len(occurrences)} groups", flush=True)
    return occurrences, molecules


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, default=DATA / "chembl37.smi")
    p.add_argument("--outdir", type=Path, default=Path("results"))
    p.add_argument("--top", type=int, default=100)
    p.add_argument(
        "--include-ring-groups",
        action="store_true",
        help="also count functional groups containing ring atoms (default: excluded)",
    )
    p.add_argument("--workers", type=int, default=None)
    a = p.parse_args()

    occ, mols = count_functional_groups(a.input, a.workers, not a.include_ring_groups)
    a.outdir.mkdir(parents=True, exist_ok=True)
    stem = a.input.stem + ("_with_ring" if a.include_ring_groups else "")
    write_counts(a.outdir / f"{stem}_functional_groups.tsv", occ, mols)
    plot_distribution(occ, mols, a.outdir / f"{stem}_fg_frequency.png", a.top)
    draw_top_structures(occ, mols, a.outdir / f"{stem}_top{a.top}_fgs.png", a.top)
    write_top_smiles(a.outdir / f"{stem}_top{a.top}_fgs.smi", occ, mols, a.top)


if __name__ == "__main__":
    main()
