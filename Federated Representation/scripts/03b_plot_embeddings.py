# 03b_plot_embeddings.py
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import torch
import umap
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


@dataclass
class Config:
    embedding_dir: Path = Path("./results/embeddings")
    plot_dir: Path = Path("./results/plots")
    meta_parquet: Path = Path("./data/processed/pgs_metadata.parquet")

    pgs_embeddings_pt: Path = Path("./results/embeddings/pgs_embeddings.pt")
    pgs_input_cosine_csv: Path = Path("./results/embeddings/pgs_input_cosine_distance.csv")

    seed: int = 42

    # UMAP params (small-n safe)
    umap_n_neighbors: int = 3
    umap_min_dist: float = 0.2
    umap_metric: str = "cosine"


def ensure_dirs(cfg: Config):
    cfg.plot_dir.mkdir(parents=True, exist_ok=True)


def pca2(X: np.ndarray) -> np.ndarray:
    Xc = X - X.mean(axis=0, keepdims=True)
    U, S, _ = np.linalg.svd(Xc, full_matrices=False)
    return U[:, :2] * S[:2]


def build_color_map(values: List[str]):
    cats = sorted(set(values))
    cmap = plt.get_cmap("tab10")
    color_map = {c: cmap(i % 10) for i, c in enumerate(cats)}
    return color_map


def plot_scatter_labeled(
    points: np.ndarray,
    labels: List[str],
    colors: List[str],
    title: str,
    outpath: Path,
):
    plt.figure(figsize=(9, 7))

    color_map = build_color_map(colors)
    for i in range(points.shape[0]):
        plt.scatter(
            points[i, 0],
            points[i, 1],
            color=color_map[colors[i]],
            s=120,
            edgecolor="black",
        )
        plt.text(
            points[i, 0],
            points[i, 1],
            labels[i],
            fontsize=9,
            ha="center",
            va="center",
        )

    legend_elements = [
        Line2D(
            [0], [0],
            marker="o",
            color="w",
            label=k,
            markerfacecolor=v,
            markersize=10,
            markeredgecolor="black",
        )
        for k, v in color_map.items()
    ]

    plt.legend(
        handles=legend_elements,
        title="Ancestry",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
    )

    plt.title(title)
    plt.xlabel("Dim 1")
    plt.ylabel("Dim 2")
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


def main():
    cfg = Config()
    ensure_dirs(cfg)

    # Load embeddings
    blob = torch.load(cfg.pgs_embeddings_pt, map_location="cpu")
    Z = blob["embeddings"].numpy()
    pgs_ids = blob["pgs_ids"]

    # Load metadata aligned to pgs_ids
    meta = pd.read_parquet(cfg.meta_parquet).copy()
    meta["pgs_id"] = meta["pgs_id"].astype(str)
    meta = meta.set_index("pgs_id").reindex(pgs_ids).reset_index()

    # Build display labels
    ancestry = meta["ancestry"].fillna("Unknown").astype(str).tolist()
    trait = meta["trait_mapped"].fillna("Unknown").astype(str).tolist()

    display_labels = [
        f"{a}\n{t}"
        for a, t in zip(ancestry, trait)
    ]

    # ---------------- Input-space cosine heatmap ----------------
    if cfg.pgs_input_cosine_csv.exists():
        dist = pd.read_csv(cfg.pgs_input_cosine_csv, index_col=0)
        plt.figure(figsize=(7, 6))
        plt.imshow(dist.to_numpy())
        plt.xticks(range(len(display_labels)), display_labels, rotation=45, ha="right")
        plt.yticks(range(len(display_labels)), display_labels)
        plt.title("PGS cosine distance (input space)")
        plt.colorbar()
        plt.tight_layout()
        plt.savefig(cfg.plot_dir / "pgs_input_cosine_distance.png")
        plt.close()

    # ---------------- PCA ----------------
    pca_xy = pca2(Z)
    plot_scatter_labeled(
        pca_xy,
        labels=display_labels,
        colors=ancestry,
        title="PGS embeddings (PCA)\nColored by ancestry, labeled by ancestry + disease",
        outpath=cfg.plot_dir / "pgs_embeddings_pca_labeled.png",
    )

    # ---------------- UMAP ----------------
    n_neighbors = min(cfg.umap_n_neighbors, max(2, len(Z) - 1))
    reducer = umap.UMAP(
        n_components=2,
        random_state=cfg.seed,
        n_neighbors=n_neighbors,
        min_dist=cfg.umap_min_dist,
        metric=cfg.umap_metric,
    )
    umap_xy = reducer.fit_transform(Z)

    plot_scatter_labeled(
        umap_xy,
        labels=display_labels,
        colors=ancestry,
        title="PGS embeddings (UMAP)\nColored by ancestry, labeled by ancestry + disease",
        outpath=cfg.plot_dir / "pgs_embeddings_umap_labeled.png",
    )

    print("Plots saved to:", cfg.plot_dir.resolve())


if __name__ == "__main__":
    main()