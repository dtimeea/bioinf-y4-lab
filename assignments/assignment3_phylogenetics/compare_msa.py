from pathlib import Path
from Bio import AlignIO, Phylo, SeqIO
from Bio.Phylo.TreeConstruction import (
    DistanceCalculator,
    DistanceTreeConstructor
)
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

folder = Path(__file__).resolve().parent

mafft = AlignIO.read(
    folder / "p53_aligned_dtimeea.fasta", "fasta"
)

clustal = AlignIO.read(
    folder / "clustalo_dtimeea.aln", "clustal"
)

# Verificarea secventelor originale
original = {
    r.id: str(r.seq).upper()
    for r in SeqIO.parse(
        folder / "p53_multi_dtimeea.fasta", "fasta"
    )
}

for alignment in [mafft, clustal]:
    sequences = {
        r.id: str(r.seq).replace("-", "").upper()
        for r in alignment
    }
    assert sequences == original, "Secventele nu coincid!"

print("Secventele originale coincid exact.")
print("MAFFT:", len(mafft), mafft.get_alignment_length())
print("Clustal Omega:", len(clustal),
      clustal.get_alignment_length())

# Construirea arborilor prin aceeasi metoda
calculator = DistanceCalculator("identity")
constructor = DistanceTreeConstructor()

dm = calculator.get_distance(clustal)
tree_clustal = constructor.nj(dm)

tree_mafft = Phylo.read(
    folder / "tree_dtimeea.nwk", "newick"
)

Phylo.write(
    tree_clustal,
    folder / "tree_clustalo_dtimeea.nwk",
    "newick"
)

pd.DataFrame(
    [[dm[a, b] for b in dm.names] for a in dm.names],
    index=dm.names,
    columns=dm.names
).to_csv(
    folder / "distance_clustalo_dtimeea.csv",
    float_format="%.6f"
)

# Compararea topologiei arborilor ne inradacinati
def splits(tree):
    all_names = frozenset(
        t.name for t in tree.get_terminals()
    )
    result = set()

    for clade in tree.find_clades():
        a = frozenset(
            t.name for t in clade.get_terminals()
        )
        b = all_names - a

        if len(a) > 1 and len(b) > 1:
            result.add(frozenset((a, b)))

    return result

s1 = splits(tree_mafft)
s2 = splits(tree_clustal)

print("\nBipartitii MAFFT:", len(s1))
print("Bipartitii Clustal:", len(s2))
print("Bipartitii comune:", len(s1 & s2))
print("Diferente topologice:", len(s1 ^ s2))

# Identificarea pozitiilor complet conservate
consensus = ""

for i in range(clustal.get_alignment_length()):
    column = clustal[:, i].upper()

    if len(set(column)) == 1 and column[0] != "-":
        consensus += column[0]
    else:
        consensus += "."

regions = []
start = None

for i, residue in enumerate(consensus + "."):
    if residue != "." and start is None:
        start = i
    elif residue == "." and start is not None:
        if i - start >= 5:
            regions.append(
                (start + 1, i, consensus[start:i], i - start)
            )
        start = None

df = pd.DataFrame(
    regions,
    columns=["Start", "End", "Sequence", "Length"]
)

df = df.sort_values("Length", ascending=False)

df.to_csv(
    folder / "conserved_regions_dtimeea.csv",
    index=False
)

print("\nRegiuni complet conservate:")
print(df.to_string(index=False))

print("\nTotal pozitii complet conservate:",
      len(consensus) - consensus.count("."))

# Aliniere cu marcarea pozitiilor conservate
with open(
    folder / "conserved_alignment_dtimeea.txt", "w"
) as f:
    for start in range(0, len(consensus), 60):
        end = min(start + 60, len(consensus))
        f.write(f"\nPositions {start+1}-{end}\n\n")

        for record in clustal:
            f.write(
                f"{record.id:22} {record.seq[start:end]}\n"
            )

        marks = "".join(
            "*" if c != "." else " "
            for c in consensus[start:end]
        )

        f.write(f"{'Conserved':22} {marks}\n")

# Vizualizarea comparativa
fig, axes = plt.subplots(1, 2, figsize=(17, 8))

for tree, ax, title in [
    (tree_mafft, axes[0], "MAFFT - Neighbor Joining"),
    (tree_clustal, axes[1], "Clustal Omega - Neighbor Joining")
]:
    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        show_confidence=False,
        label_func=lambda c: c.name if c.is_terminal() else None
    )
    ax.set_title(title)

plt.tight_layout()
fig.savefig(
    folder / "tree_comparison_dtimeea.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close(fig)

print("\nComparatia a fost salvata.")
