"""Count element frequencies (atoms per element) over a SMILES file."""

import argparse
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger

from moltally.common import (
    CHUNK,
    DATA,
    PROGRESS_EVERY,
    parallel_map,
    read_chunks,
)


def _count_chunk(lines: list[str]) -> tuple[Counter, Counter, int]:
    RDLogger.DisableLog("rdApp.*")
    atoms, molecules, n_parsed = Counter(), Counter(), 0
    for line in lines:
        mol = Chem.MolFromSmiles(line.split()[0])
        if mol is None:
            continue
        n_parsed += 1
        symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
        atoms.update(symbols)
        molecules.update(set(symbols))
    return atoms, molecules, n_parsed


def count_elements(
    path: Path, workers: int | None = None
) -> tuple[Counter, Counter, int]:
    """Atoms per element, molecules per element, and the number of parsed molecules."""
    atoms, molecules, total = Counter(), Counter(), 0
    for i, (a, m, n) in enumerate(parallel_map(_count_chunk, read_chunks(path), workers), 1):
        atoms.update(a)
        molecules.update(m)
        total += n
        if i % PROGRESS_EVERY == 0:
            print(f"{i * CHUNK:>11,} molecules, {len(atoms)} elements", flush=True)
    return atoms, molecules, total


def ranked_elements(atoms: Counter, molecules: Counter) -> list[tuple[str, int, int]]:
    """(element, molecules, atoms) sorted by molecules, then atoms, then symbol."""
    items = [(el, molecules[el], n) for el, n in atoms.items()]
    return sorted(items, key=lambda t: (-t[1], -t[2], t[0]))


def write_table(path: Path, rows: list[tuple[str, int, int]], n_molecules: int) -> None:
    """pct_molecules: share of molecules containing the element;
    pct_atoms: share of all atoms that are this element."""
    total_atoms = sum(r[2] for r in rows)
    with open(path, "w") as fh:
        fh.write("element\tmolecules\tpct_molecules\tatoms\tpct_atoms\n")
        for el, n_mol, n_atoms in rows:
            fh.write(
                f"{el}\t{n_mol}\t{100 * n_mol / n_molecules:.4f}\t"
                f"{n_atoms}\t{100 * n_atoms / total_atoms:.4f}\n"
            )


def plot_elements(rows: list[tuple[str, int, int]], path: Path) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(max(6, 0.5 * len(rows) + 2), 4.5))
    x = range(len(rows))
    ax.bar([i - 0.2 for i in x], [r[1] for r in rows], width=0.4, label="molecules")
    ax.bar([i + 0.2 for i in x], [r[2] for r in rows], width=0.4, label="atoms")
    ax.set_xticks(list(x), [r[0] for r in rows])
    ax.set_yscale("log")
    ax.set(xlabel="element", ylabel="count", title=f"Element frequencies ({len(rows)} elements)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, default=DATA / "chembl37.smi")
    p.add_argument("--outdir", type=Path, default=Path("results"))
    p.add_argument("--workers", type=int, default=None)
    a = p.parse_args()

    atoms, molecules, total = count_elements(a.input, a.workers)
    rows = ranked_elements(atoms, molecules)
    a.outdir.mkdir(parents=True, exist_ok=True)
    stem = a.input.stem
    write_table(a.outdir / f"{stem}_elements.tsv", rows, total)
    plot_elements(rows, a.outdir / f"{stem}_element_frequency.png")


if __name__ == "__main__":
    main()
