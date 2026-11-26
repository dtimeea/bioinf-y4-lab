from pathlib import Path
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import (
    DistanceCalculator,
    DistanceTreeConstructor
)
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

folder = Path(__file__).resolve().parent

# 1. Citirea alinierii multiple
alignment = AlignIO.read(
    folder / "p53_aligned_dtimeea.fasta",
    "fasta"
)

assert len(alignment) >= 10

# 2. Matricea de distante
calculator = DistanceCalculator("identity")
dm = calculator.get_distance(alignment)

df = pd.DataFrame(
    [[dm[a, b] for b in dm.names] for a in dm.names],
    index=dm.names,
    columns=dm.names
)

df.to_csv(
    folder / "distance_matrix_dtimeea.csv",
    float_format="%.6f"
)

print("Matricea de distante:")
print(df.round(4).to_string())

# 3. Perechile extreme
pairs = [
    (dm[a, b], a, b)
    for i, a in enumerate(dm.names)
    for b in dm.names[i + 1:]
]

closest = min(pairs)
farthest = max(pairs)

print("\nPerechea cea mai apropiata:", closest)
print("Perechea cea mai indepartata:", farthest)

# 4. Arbore Neighbor-Joining
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)

Phylo.write(
    tree,
    folder / "tree_dtimeea.nwk",
    "newick"
)

print("\nArbore Neighbor-Joining:")
Phylo.draw_ascii(tree)

# 5. Vizualizare grafica - bonus
fig, ax = plt.subplots(figsize=(13, 8))

Phylo.draw(
    tree,
    axes=ax,
    do_show=False,
    show_confidence=False,
    label_func=lambda c: c.name if c.is_terminal() else None
)

ax.set_title(
    "p53 protein phylogeny - Neighbor-Joining\n"
    "11 species | MAFFT alignment"
)

ax.set_xlabel("Branch length (identity distance)")
ax.set_ylabel("")

plt.tight_layout()

fig.savefig(
    folder / "tree_dtimeea.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)

print("\nMatricea, arborele si imaginea au fost salvate.")
