from __future__ import annotations
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
import networkx as nx


# Configurare
INPUT_CSV = Path("data/work/dtimeea/lab06/expression_matrix.csv")
OUTPUT_DIR = Path("labs/06_wgcna/submissions/dtimeea")
OUTPUT_CSV = OUTPUT_DIR / "modules_dtimeea.csv"

CORR_METHOD = "spearman"
VARIANCE_THRESHOLD = 0.1
ADJ_THRESHOLD = 0.7
USE_ABS_CORR = False
MAKE_UNDIRECTED = True
ALREADY_LOG_TRANSFORMED = True


def read_expression_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna()
    df = df.groupby(level=0).mean()
    return df


def log_and_filter(df: pd.DataFrame,
                   variance_threshold: float) -> pd.DataFrame:
    # GSE96058 contine valori deja transformate logaritmic.
    if not ALREADY_LOG_TRANSFORMED:
        if (df < 0).any().any():
            raise ValueError("Datele brute trebuie sa fie nenegative.")
        df = np.log2(df + 1)

    variances = df.var(axis=1)
    return df.loc[variances >= variance_threshold]


def correlation_matrix(df: pd.DataFrame,
                       method: str = "spearman",
                       use_abs: bool = False) -> pd.DataFrame:
    corr = df.T.corr(method=method)
    corr = corr.fillna(0)

    if use_abs:
        corr = corr.abs()

    return corr


def adjacency_from_correlation(corr: pd.DataFrame,
                               threshold: float,
                               weighted: bool = False) -> pd.DataFrame:
    if weighted:
        A = corr.where(corr >= threshold, 0)
    else:
        A = (corr >= threshold).astype(int)

    np.fill_diagonal(A.values, 0)
    return A


def graph_from_adjacency(A: pd.DataFrame,
                         undirected: bool = True) -> nx.Graph:
    if undirected:
        G = nx.from_pandas_adjacency(A)
    else:
        G = nx.from_pandas_adjacency(
            A, create_using=nx.DiGraph
        )

    isolates = list(nx.isolates(G))
    print("Gene izolate eliminate:", len(isolates))
    G.remove_nodes_from(isolates)

    return G


def detect_modules_louvain_or_greedy(
        G: nx.Graph) -> Dict[str, int]:

    if G.number_of_edges() == 0:
        return {}

    if hasattr(nx.community, "louvain_communities"):
        communities = nx.community.louvain_communities(
            G, seed=42
        )
        print("Algoritm utilizat: Louvain")
    else:
        communities = nx.community.greedy_modularity_communities(G)
        print("Algoritm utilizat: Greedy modularity")

    mapping = {}

    for module_id, community in enumerate(communities, start=1):
        for gene in community:
            mapping[gene] = module_id

    return mapping


def save_modules_csv(mapping: Dict[str, int],
                     out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    df_modules = pd.DataFrame({
        "Gene": list(mapping.keys()),
        "Module": list(mapping.values())
    }).sort_values(["Module", "Gene"])

    df_modules.to_csv(out_csv, index=False)


if __name__ == "__main__":
    df = read_expression_matrix(INPUT_CSV)
    print("Dimensiuni matrice initiala:", df.shape)

    df = log_and_filter(df, VARIANCE_THRESHOLD)
    print("Dimensiuni dupa filtrare:", df.shape)

    corr = correlation_matrix(
        df, CORR_METHOD, USE_ABS_CORR
    )

    A = adjacency_from_correlation(
        corr, ADJ_THRESHOLD
    )

    G = graph_from_adjacency(A, MAKE_UNDIRECTED)

    print("Metrica de corelatie:", CORR_METHOD)
    print("Prag de corelatie:", ADJ_THRESHOLD)
    print("Numar noduri:", G.number_of_nodes())
    print("Numar muchii:", G.number_of_edges())

    gene_to_module = detect_modules_louvain_or_greedy(G)

    print("Numar module:", len(set(gene_to_module.values())))

    save_modules_csv(gene_to_module, OUTPUT_CSV)

    print("\nDistributia genelor pe module:")
    print(pd.Series(
        list(gene_to_module.values()), dtype=int
    ).value_counts().sort_index())

    print("\nFisier salvat:", OUTPUT_CSV)
