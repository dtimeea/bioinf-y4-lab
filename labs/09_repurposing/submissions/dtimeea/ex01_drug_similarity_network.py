from pathlib import Path
from itertools import combinations

import networkx as nx
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HANDLE = "dtimeea"

DRUG_GENE_CSV = Path(
    f"data/work/{HANDLE}/lab09/drug_gene_{HANDLE}.csv"
)

OUT_DIR = Path(
    f"labs/09_repurposing/submissions/{HANDLE}"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_DRUG_SUMMARY = OUT_DIR / f"drug_summary_{HANDLE}.csv"
OUT_DRUG_SIMILARITY = OUT_DIR / f"drug_similarity_{HANDLE}.csv"


def ensure_exists(path):
    if not path.is_file():
        raise FileNotFoundError(f"Missing file: {path}")


def load_drug_gene_table(path):
    ensure_exists(path)
    df = pd.read_csv(path, comment="#")

    if not {"drug", "gene"}.issubset(df.columns):
        raise ValueError("Required columns: drug, gene")

    df = df.dropna(subset=["drug", "gene"]).copy()
    df["drug"] = df["drug"].str.strip()
    df["gene"] = df["gene"].str.strip()

    return df[
        (df["drug"] != "") & (df["gene"] != "")
    ].drop_duplicates(["drug", "gene"])


def build_drug2genes(df):
    return df.groupby("drug")["gene"].apply(set).to_dict()


def build_bipartite_graph(drug2genes):
    B = nx.Graph()

    for drug, genes in sorted(drug2genes.items()):
        B.add_node(drug, bipartite="drug")

        for gene in sorted(genes):
            B.add_node(gene, bipartite="gene")
            B.add_edge(drug, gene)

    return B


def summarize_drugs(drug2genes):
    rows = [
        {"drug": drug, "num_targets": len(genes)}
        for drug, genes in drug2genes.items()
    ]

    return pd.DataFrame(rows).sort_values(
        ["num_targets", "drug"],
        ascending=[False, True]
    )


def jaccard_similarity(s1, s2):
    union = s1 | s2
    return len(s1 & s2) / len(union) if union else 0.0


def compute_drug_similarity_edges(drug2genes, min_sim=0.0):
    edges = []

    for d1, d2 in combinations(sorted(drug2genes), 2):
        similarity = jaccard_similarity(
            drug2genes[d1], drug2genes[d2]
        )

        if similarity >= min_sim:
            edges.append((d1, d2, similarity))

    return edges


def edges_to_dataframe(edges):
    return pd.DataFrame(
        edges,
        columns=["drug1", "drug2", "similarity"]
    ).sort_values(
        ["similarity", "drug1", "drug2"],
        ascending=[False, True, True]
    )


def plot_bipartite(B, summary):
    top_drugs = summary.head(20)["drug"].tolist()

    genes = [
        n for n, d in B.nodes(data=True)
        if d["bipartite"] == "gene"
    ]

    H = B.subgraph(top_drugs + genes).copy()

    pos = nx.bipartite_layout(H, top_drugs)

    colors = [
        "royalblue" if H.nodes[n]["bipartite"] == "drug"
        else "tomato"
        for n in H.nodes()
    ]

    sizes = [300 + 120 * H.degree(n) for n in H.nodes()]

    plt.figure(figsize=(16, 13))

    nx.draw(
        H, pos,
        node_color=colors,
        node_size=sizes,
        with_labels=True,
        font_size=8,
        edge_color="gray",
        alpha=0.85
    )

    plt.title(
        "Drug-Gene Network - DGIdb 2026-10\n"
        "Top 20 drugs by number of associated genes"
    )

    plt.tight_layout()
    plt.savefig(
        OUT_DIR / f"network_drug_gene_{HANDLE}.png",
        dpi=300, bbox_inches="tight"
    )
    plt.close()


def plot_similarity(G, summary):
    top_drugs = summary.head(20)["drug"].tolist()
    H = G.subgraph(top_drugs).copy()

    pos = nx.spring_layout(H, seed=42, weight="weight")

    widths = [
        0.5 + 4 * H[u][v]["weight"]
        for u, v in H.edges()
    ]

    plt.figure(figsize=(16, 12))

    nx.draw(
        H, pos,
        with_labels=True,
        node_color="skyblue",
        node_size=1000,
        font_size=8,
        width=widths,
        edge_color="gray",
        alpha=0.8
    )

    plt.title(
        "Drug Similarity Network - Jaccard\n"
        "Top 20 drugs by number of associated genes"
    )

    plt.tight_layout()
    plt.savefig(
        OUT_DIR / f"drug_similarity_{HANDLE}.png",
        dpi=300, bbox_inches="tight"
    )
    plt.close()


if __name__ == "__main__":
    df = load_drug_gene_table(DRUG_GENE_CSV)

    drug2genes = build_drug2genes(df)
    B = build_bipartite_graph(drug2genes)

    summary = summarize_drugs(drug2genes)
    summary.to_csv(OUT_DRUG_SUMMARY, index=False)

    edges = compute_drug_similarity_edges(drug2genes)
    similarity = edges_to_dataframe(edges)

    similarity.to_csv(OUT_DRUG_SIMILARITY, index=False)

    similarity.to_csv(
        OUT_DIR / f"REPURPOSING_{HANDLE}.csv",
        index=False
    )

    G = nx.Graph()
    G.add_nodes_from(drug2genes)

    G.add_weighted_edges_from(
        (d1, d2, w)
        for d1, d2, w in edges
        if w > 0
    )

    plot_bipartite(B, summary)
    plot_similarity(G, summary)

    print("Medicamente:", len(drug2genes))
    print("Gene:", df["gene"].nunique())
    print("Noduri bipartite:", B.number_of_nodes())
    print("Muchii bipartite:", B.number_of_edges())
    print("Componente conexe:", nx.number_connected_components(B))
    print("Perechi evaluate:", len(similarity))
    print("Muchii similarity:", G.number_of_edges())

    print("\nTop 10 medicamente:")
    print(summary.head(10).to_string(index=False))

    print("\nTop 10 similaritati:")
    print(similarity.head(10).to_string(index=False))

    print("\nExercise 9.1 completed.")
