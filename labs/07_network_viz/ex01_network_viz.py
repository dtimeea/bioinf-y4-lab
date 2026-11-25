from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# Configurare
HANDLE = "dtimeea"

EXPR_CSV = Path(
    f"data/work/{HANDLE}/lab06/expression_matrix.csv"
)

MODULES_CSV = Path(
    f"labs/06_wgcna/submissions/{HANDLE}/modules_{HANDLE}.csv"
)

PRECOMPUTED_ADJ_CSV: Optional[Path] = None

CORR_METHOD = "spearman"
USE_ABS_CORR = False
ADJ_THRESHOLD = 0.7
VARIANCE_THRESHOLD = 0.1
WEIGHTED = False

SEED = 42
TOPK_HUBS = 10
NODE_BASE_SIZE = 300
EDGE_ALPHA = 0.4

OUT_DIR = Path(
    f"labs/07_network_viz/submissions/{HANDLE}"
)

OUT_PNG = OUT_DIR / f"network_{HANDLE}.png"
OUT_HUBS = OUT_DIR / f"hubs_{HANDLE}.csv"


def ensure_exists(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)


def read_expression_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna()
    df = df.groupby(level=0).mean()
    return df.loc[df.var(axis=1) >= VARIANCE_THRESHOLD]


def read_modules_csv(path: Path) -> Dict[str, int]:
    df = pd.read_csv(path)
    return dict(zip(df["Gene"], df["Module"]))


def correlation_to_adjacency(
        expr: pd.DataFrame,
        method: str,
        use_abs: bool,
        threshold: float,
        weighted: bool) -> pd.DataFrame:

    corr = expr.T.corr(method=method).fillna(0)

    if use_abs:
        corr = corr.abs()

    if weighted:
        A = corr.where(corr >= threshold, 0)
    else:
        A = (corr >= threshold).astype(int)

    np.fill_diagonal(A.values, 0)

    return A


def graph_from_adjacency(A: pd.DataFrame) -> nx.Graph:
    G = nx.from_pandas_adjacency(A)
    G.remove_nodes_from(list(nx.isolates(G)))
    return G


def color_map_from_modules(
        nodes: Iterable[str],
        gene2module: Dict[str, int]) -> dict:

    modules = sorted(set(gene2module.values()))
    palette = plt.get_cmap("tab10")

    module_colors = {
        module: palette(i % 10)
        for i, module in enumerate(modules)
    }

    return {
        gene: module_colors[gene2module[gene]]
        for gene in nodes
    }


def compute_hubs(G: nx.Graph, topk: int) -> pd.DataFrame:
    degrees = dict(G.degree())
    betweenness = nx.betweenness_centrality(G)

    df = pd.DataFrame({
        "Gene": list(G.nodes()),
        "Degree": [degrees[g] for g in G.nodes()],
        "Betweenness": [
            betweenness[g] for g in G.nodes()
        ]
    })

    return df.sort_values(
        ["Degree", "Betweenness", "Gene"],
        ascending=[False, False, True]
    ).head(topk)


if __name__ == "__main__":

    ensure_exists(EXPR_CSV)
    ensure_exists(MODULES_CSV)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    expr = read_expression_matrix(EXPR_CSV)
    gene2module = read_modules_csv(MODULES_CSV)

    if PRECOMPUTED_ADJ_CSV is not None:
        ensure_exists(PRECOMPUTED_ADJ_CSV)
        A = pd.read_csv(
            PRECOMPUTED_ADJ_CSV, index_col=0
        )
    else:
        A = correlation_to_adjacency(
            expr,
            CORR_METHOD,
            USE_ABS_CORR,
            ADJ_THRESHOLD,
            WEIGHTED
        )

    G = graph_from_adjacency(A)

    # Verificam concordanta cu modulele din Lab 6
    if set(G.nodes()) != set(gene2module):
        raise ValueError(
            "Nodurile nu coincid cu modulele din Lab 6."
        )

    print("Numar noduri:", G.number_of_nodes())
    print("Numar muchii:", G.number_of_edges())
    print("Numar module:", len(set(gene2module.values())))

    colors = color_map_from_modules(
        G.nodes(), gene2module
    )

    hubs_df = compute_hubs(G, TOPK_HUBS)
    hubs_df["Module"] = hubs_df["Gene"].map(gene2module)

    hubs = set(hubs_df["Gene"])

    node_sizes = [
        NODE_BASE_SIZE + 100 * G.degree(g)
        for g in G.nodes()
    ]

    node_colors = [
        colors[g] for g in G.nodes()
    ]

    borders = [
        "black" if g in hubs else "white"
        for g in G.nodes()
    ]

    pos = nx.spring_layout(
        G, seed=SEED, k=0.9, iterations=200
    )

    fig, ax = plt.subplots(figsize=(12, 9))

    nx.draw_networkx_edges(
        G, pos,
        alpha=EDGE_ALPHA,
        edge_color="gray",
        ax=ax
    )

    nx.draw_networkx_nodes(
        G, pos,
        node_color=node_colors,
        node_size=node_sizes,
        edgecolors=borders,
        linewidths=1.8,
        ax=ax
    )

    nx.draw_networkx_labels(
        G, pos,
        labels={g: g for g in hubs},
        font_size=8,
        ax=ax
    )

    legend = []

    for module in sorted(set(gene2module.values())):
        gene = next(
            g for g in G.nodes()
            if gene2module[g] == module
        )

        legend.append(Line2D(
            [0], [0],
            marker="o",
            color="w",
            label=f"Module {module}",
            markerfacecolor=colors[gene],
            markersize=10
        ))

    legend.append(Line2D(
        [0], [0],
        marker="o",
        color="w",
        label="Top 10 hub genes",
        markerfacecolor="white",
        markeredgecolor="black",
        markersize=10
    ))

    ax.legend(handles=legend, loc="upper right")

    ax.set_title(
        "Gene Co-expression Network - GSE96058\n"
        "Spearman >= 0.7 | Louvain modules",
        fontsize=14
    )

    ax.axis("off")
    fig.tight_layout()

    fig.savefig(
        OUT_PNG, dpi=300, bbox_inches="tight"
    )

    plt.close(fig)

    hubs_df.to_csv(OUT_HUBS, index=False)

    print("\nTop hub genes:")
    print(hubs_df.to_string(index=False))

    print("\nImagine salvata:", OUT_PNG)
    print("Hub genes salvate:", OUT_HUBS)
