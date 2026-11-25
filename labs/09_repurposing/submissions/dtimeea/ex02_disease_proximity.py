from pathlib import Path

import networkx as nx
import pandas as pd

from ex01_drug_similarity_network import (
    load_drug_gene_table,
    build_drug2genes,
    build_bipartite_graph
)

HANDLE = "dtimeea"

DRUG_GENE_CSV = Path(
    f"data/work/{HANDLE}/lab09/drug_gene_{HANDLE}.csv"
)

DISEASE_GENES_TXT = Path(
    f"data/work/{HANDLE}/lab09/disease_genes_{HANDLE}.txt"
)

OUT_DIR = Path(
    f"labs/09_repurposing/submissions/{HANDLE}"
)

OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_DRUG_PRIORITY = (
    OUT_DIR / f"drug_priority_{HANDLE}.csv"
)


def ensure_exists(path):
    if not path.is_file():
        raise FileNotFoundError(f"Missing file: {path}")


def load_bipartite_graph_or_build():
    df = load_drug_gene_table(DRUG_GENE_CSV)
    drug2genes = build_drug2genes(df)

    return build_bipartite_graph(drug2genes)


def load_disease_genes(path):
    ensure_exists(path)

    return {
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    }


def get_drug_nodes(B):
    return sorted(
        n for n, d in B.nodes(data=True)
        if d.get("bipartite") == "drug"
    )


def compute_drug_disease_distance(
    B, drug, disease_genes,
    mode="mean", max_dist=5
):
    if mode not in {"mean", "min"}:
        raise ValueError("Invalid distance mode")

    lengths = nx.single_source_shortest_path_length(
        B, drug
    )

    distances = []

    for gene in sorted(disease_genes):
        if gene not in B:
            continue

        distances.append(
            lengths.get(
                gene,
                max(max_dist + 1, B.number_of_nodes())
            )
        )

    if not distances:
        raise ValueError("No disease genes in graph")

    if mode == "min":
        return float(min(distances))

    return sum(distances) / len(distances)


def rank_drugs_by_proximity(
    B, disease_genes, mode="mean"
):
    rows = []

    for drug in get_drug_nodes(B):
        distance = compute_drug_disease_distance(
            B, drug, disease_genes, mode
        )

        rows.append({
            "drug": drug,
            "distance": distance
        })

    return pd.DataFrame(rows).sort_values(
        ["distance", "drug"],
        ascending=[True, True]
    ).reset_index(drop=True)


if __name__ == "__main__":
    ensure_exists(DRUG_GENE_CSV)
    ensure_exists(DISEASE_GENES_TXT)

    B = load_bipartite_graph_or_build()

    disease_genes = load_disease_genes(
        DISEASE_GENES_TXT
    )

    present = {
        gene for gene in disease_genes
        if gene in B
        and B.nodes[gene].get("bipartite") == "gene"
    }

    missing = disease_genes - present

    print("Disease genes:", sorted(disease_genes))
    print("Gene prezente:", sorted(present))
    print("Gene absente:", sorted(missing))

    if not present:
        raise ValueError("No disease genes found")

    ranking = rank_drugs_by_proximity(
        B, present, mode="mean"
    )

    ranking.to_csv(
        OUT_DRUG_PRIORITY, index=False
    )

    print("\nMedicamente evaluate:", len(ranking))

    print("\nTop 15 medicamente dupa proximitate:")
    print(ranking.head(15).to_string(index=False))

    print("\nDistanta minima:", ranking["distance"].min())
    print("Distanta maxima:", ranking["distance"].max())
    print("Distanta medie:", ranking["distance"].mean())

    print("\nFisier salvat:", OUT_DRUG_PRIORITY)
    print("Exercise 9.2 completed.")
