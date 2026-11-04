"""Shared helpers: chunked reading, count tables, frequency plot, structure grid."""

from collections import Counter, deque
from itertools import islice
from multiprocessing import Pool
from pathlib import Path

from rdkit import Chem

DATA = Path(__file__).parent / "data"
CHUNK = 5000
PROGRESS_EVERY = 200  # chunks between progress messages (1M molecules)


def ranked(
    occurrences: Counter, molecules: Counter, top: int | None = None
) -> list[tuple[str, int, int]]:
    """(smiles, molecules, occurrences) sorted by number of molecules, descending.

    Ties are broken by occurrences, then alphabetically, for reproducibility.
    """
    items = [(smi, molecules[smi], n) for smi, n in occurrences.items()]
    items.sort(key=lambda t: (-t[1], -t[2], t[0]))
    return items[:top]


def read_chunks(path: Path):
    """Yield lists of up to CHUNK non-empty lines of a SMILES file."""
    with open(path) as fh:
        lines = (line for line in fh if line.strip())
        while chunk := list(islice(lines, CHUNK)):
            yield chunk


def parallel_map(func, chunks, workers: int | None = None, prefetch: int = 4):
    """Ordered parallel map over a (possibly huge) iterable with bounded memory.

    Unlike Pool.imap, which eagerly reads the whole input iterable into its task
    queue, at most ``prefetch * workers`` chunks are in flight at any time.
    """
    with Pool(workers) as pool:
        limit = prefetch * pool._processes
        pending: deque = deque()
        for chunk in chunks:
            pending.append(pool.apply_async(func, (chunk,)))
            if len(pending) >= limit:
                yield pending.popleft().get()
        while pending:
            yield pending.popleft().get()


def write_counts(path: Path, occurrences: Counter, molecules: Counter) -> None:
    with open(path, "w") as fh:
        fh.write("smiles\tmolecules\toccurrences\n")
        for smi, n_mol, n_occ in ranked(occurrences, molecules):
            fh.write(f"{smi}\t{n_mol}\t{n_occ}\n")


def write_top_smiles(path: Path, occurrences: Counter, molecules: Counter, top: int) -> None:
    path.write_text("".join(f"{smi}\n" for smi, _, _ in ranked(occurrences, molecules, top)))


def plot_distribution(
    occurrences: Counter,
    molecules: Counter,
    path: Path,
    top: int = 100,
    label: str = "functional groups",
) -> None:
    import matplotlib.pyplot as plt

    counts = [n for _, n, _ in ranked(occurrences, molecules)]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.5))
    ax1.loglog(range(1, len(counts) + 1), counts, ".", ms=3)
    ax1.axvline(top, color="gray", ls="--", lw=1, label=f"top {top}")
    ax1.set(xlabel="rank", ylabel="molecules", title=f"Rank-frequency ({len(counts)} {label})")
    ax1.legend()
    ax2.bar(range(1, top + 1), counts[:top])
    ax2.set_yscale("log")
    ax2.set(xlabel="rank", ylabel="molecules", title=f"Top {top} {label}")
    fig.tight_layout()
    fig.savefig(path, dpi=150)


def draw_top_structures(
    occurrences: Counter, molecules: Counter, path: Path, top: int = 100, per_row: int = 10
) -> None:
    """Grid image of the most common items, labelled with rank, SMILES and molecule count."""
    import matplotlib.pyplot as plt
    from rdkit.Chem import Draw

    items = ranked(occurrences, molecules, top)
    rows = -(-len(items) // per_row)
    fig, axes = plt.subplots(rows, per_row, figsize=(per_row * 2.4, rows * 2.4), squeeze=False)
    for ax in axes.flat:
        ax.axis("off")
    for ax, (rank, (smi, n, _)) in zip(axes.flat, enumerate(items, 1)):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:  # fragments are not always valid sanitized molecules
            mol = Chem.MolFromSmiles(smi, sanitize=False)
            mol.UpdatePropertyCache(strict=False)
            Chem.FastFindRings(mol)
        ax.imshow(Draw.MolToImage(mol, size=(300, 240)))
        ax.set_title(f"#{rank} {smi}\n{n:,}", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
