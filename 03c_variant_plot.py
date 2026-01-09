# 03c_plot_variant_embeddings_by_ancestry_trait_sharedness.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import numpy as np
import pandas as pd
import umap
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


@dataclass
class Config:
    variant_embeddings_csv: Path = Path("./results/embeddings/variant_embeddings.csv")
    master_annotated_parquet: Path = Path("./data/processed/pgs_variants_master_annotated.parquet")
    plot_dir: Path = Path("./results/plots")

    seed: int = 42
    umap_n_neighbors: int = 20
    umap_min_dist: float = 0.15
    umap_metric: str = "cosine"

    label_top_n: int = 30  # label top variants by mean_abs_weight
    shared_threshold: int = 2  # shared if present in >=2 ancestries


def ensure_dirs(cfg: Config):
    cfg.plot_dir.mkdir(parents=True, exist_ok=True)


def build_color_map(values, palette="tab20"):
    cats = sorted(pd.Series(values).fillna("Unknown").astype(str).unique().tolist())
    cmap = plt.get_cmap(palette)
    return {c: cmap(i % cmap.N) for i, c in enumerate(cats)}


def scatter(points, categories, title, outpath, legend_title="Category", labels=None):
    categories = pd.Series(categories).fillna("Unknown").astype(str).tolist()
    cmap = build_color_map(categories, "tab20")

    plt.figure(figsize=(10, 8))
    for i in range(points.shape[0]):
        plt.scatter(points[i, 0], points[i, 1], s=18, color=cmap[categories[i]], alpha=0.85)
        if labels is not None and labels[i]:
            plt.text(points[i, 0], points[i, 1], labels[i], fontsize=7)

    legend_elems = [
        Line2D([0], [0], marker="o", color="w", label=k, markerfacecolor=v, markersize=8)
        for k, v in cmap.items()
    ]
    plt.legend(handles=legend_elems, title=legend_title, bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.title(title)
    plt.xlabel("Dim 1")
    plt.ylabel("Dim 2")
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


def pca2(X: np.ndarray) -> np.ndarray:
    Xc = X - X.mean(axis=0, keepdims=True)
    U, S, _ = np.linalg.svd(Xc, full_matrices=False)
    return U[:, :2] * S[:2]


def main():
    cfg = Config()
    ensure_dirs(cfg)

    vdf = pd.read_csv(cfg.variant_embeddings_csv)
    z_cols = [c for c in vdf.columns if c.startswith("z")]
    if len(z_cols) < 2:
        raise ValueError("Need at least 2 latent dims in variant_embeddings.csv")

    # Load master to compute ancestry/trait presence per locus_id
    m = pd.read_parquet(cfg.master_annotated_parquet).copy()
    m["locus_id"] = m["locus_id"].astype(str)

    if "weight_scaled" not in m.columns:
        raise ValueError("master_annotated parquet must contain weight_scaled")

    m["absw"] = pd.to_numeric(m["weight_scaled"], errors="coerce").fillna(0).abs()
    m["ancestry"] = m.get("ancestry", "Unknown").fillna("Unknown").astype(str)
    m["trait_mapped"] = m.get("trait_mapped", "Unknown").fillna("Unknown").astype(str)

    # ---- primary ancestry per locus (largest sum abs weight) ----
    anc_score = m.groupby(["locus_id", "ancestry"])["absw"].sum().reset_index()
    anc_best = anc_score.sort_values(["locus_id", "absw"], ascending=[True, False]).drop_duplicates("locus_id")
    anc_best = anc_best.rename(columns={"ancestry": "primary_ancestry", "absw": "anc_absw_sum"})[
        ["locus_id", "primary_ancestry", "anc_absw_sum"]
    ]

    # ---- primary trait per locus ----
    tr_score = m.groupby(["locus_id", "trait_mapped"])["absw"].sum().reset_index()
    tr_best = tr_score.sort_values(["locus_id", "absw"], ascending=[True, False]).drop_duplicates("locus_id")
    tr_best = tr_best.rename(columns={"trait_mapped": "primary_trait", "absw": "trait_absw_sum"})[
        ["locus_id", "primary_trait", "trait_absw_sum"]
    ]

    # ---- sharedness across ancestries ----
    anc_presence = (
        m.loc[m["absw"] > 0, ["locus_id", "ancestry"]]
        .drop_duplicates()
        .groupby("locus_id")["ancestry"]
        .nunique()
        .reset_index()
        .rename(columns={"ancestry": "n_ancestries"})
    )

    def shared_label(n):
        n = int(n)
        if n <= 1:
            return "ancestry_specific"
        if n == 2:
            return "shared_2"
        if n == 3:
            return "shared_3"
        if n == 4:
            return "shared_4"
        return "shared_5plus"

    anc_presence["sharedness"] = anc_presence["n_ancestries"].apply(shared_label)
    anc_presence["shared_binary"] = anc_presence["n_ancestries"].apply(
        lambda x: "shared_2plus" if int(x) >= cfg.shared_threshold else "ancestry_specific"
    )

    # Merge all locus-level summaries onto variant embedding table
    vdf["locus_id"] = vdf["locus_id"].astype(str)
    vdf = (
        vdf.merge(anc_best, on="locus_id", how="left")
           .merge(tr_best, on="locus_id", how="left")
           .merge(anc_presence, on="locus_id", how="left")
    )

    vdf["primary_ancestry"] = vdf["primary_ancestry"].fillna("Unknown")
    vdf["primary_trait"] = vdf["primary_trait"].fillna("Unknown")
    vdf["n_ancestries"] = pd.to_numeric(vdf["n_ancestries"], errors="coerce").fillna(0).astype(int)
    vdf["sharedness"] = vdf["sharedness"].fillna("Unknown")
    vdf["shared_binary"] = vdf["shared_binary"].fillna("Unknown")

    # Labels for top variants
    vdf["mean_abs_weight"] = pd.to_numeric(vdf.get("mean_abs_weight", 0), errors="coerce").fillna(0)
    top = set(vdf.sort_values("mean_abs_weight", ascending=False).head(cfg.label_top_n)["locus_id"].tolist())
    labels = [lid if lid in top else "" for lid in vdf["locus_id"].tolist()]

    Z = vdf[z_cols].to_numpy().astype(np.float32)

    # PCA plots
    P = pca2(Z)
    scatter(P, vdf["primary_ancestry"], "Variant embeddings (PCA) colored by PRIMARY ancestry",
            cfg.plot_dir / "variant_pca_by_primary_ancestry.png",
            legend_title="Primary ancestry", labels=labels)

    scatter(P, vdf["primary_trait"], "Variant embeddings (PCA) colored by PRIMARY trait",
            cfg.plot_dir / "variant_pca_by_primary_trait.png",
            legend_title="Primary trait", labels=labels)

    scatter(P, vdf["shared_binary"], f"Variant embeddings (PCA) colored by SHARED vs SPECIFIC (>= {cfg.shared_threshold} ancestries)",
            cfg.plot_dir / "variant_pca_by_shared_binary.png",
            legend_title="Sharedness", labels=labels)

    scatter(P, vdf["sharedness"], "Variant embeddings (PCA) colored by # ancestries bin",
            cfg.plot_dir / "variant_pca_by_sharedness_bins.png",
            legend_title="Ancestry count bin", labels=labels)

    scatter(P, vdf["n_ancestries"].astype(str), "Variant embeddings (PCA) colored by # ancestries (exact)",
            cfg.plot_dir / "variant_pca_by_n_ancestries_exact.png",
            legend_title="# ancestries", labels=labels)

    # UMAP plots
    n_neighbors = min(cfg.umap_n_neighbors, max(2, len(vdf) - 1))
    reducer = umap.UMAP(
        n_components=2,
        random_state=cfg.seed,
        n_neighbors=n_neighbors,
        min_dist=cfg.umap_min_dist,
        metric=cfg.umap_metric,
    )
    U = reducer.fit_transform(Z)

    scatter(U, vdf["primary_ancestry"], f"Variant embeddings (UMAP) colored by PRIMARY ancestry (n_neighbors={n_neighbors})",
            cfg.plot_dir / "variant_umap_by_primary_ancestry.png",
            legend_title="Primary ancestry", labels=labels)

    scatter(U, vdf["primary_trait"], f"Variant embeddings (UMAP) colored by PRIMARY trait (n_neighbors={n_neighbors})",
            cfg.plot_dir / "variant_umap_by_primary_trait.png",
            legend_title="Primary trait", labels=labels)

    scatter(U, vdf["shared_binary"], f"Variant embeddings (UMAP) colored by SHARED vs SPECIFIC (>= {cfg.shared_threshold} ancestries)",
            cfg.plot_dir / "variant_umap_by_shared_binary.png",
            legend_title="Sharedness", labels=labels)

    scatter(U, vdf["sharedness"], f"Variant embeddings (UMAP) colored by # ancestries bin (n_neighbors={n_neighbors})",
            cfg.plot_dir / "variant_umap_by_sharedness_bins.png",
            legend_title="Ancestry count bin", labels=labels)

    scatter(U, vdf["n_ancestries"].astype(str), f"Variant embeddings (UMAP) colored by # ancestries (exact) (n_neighbors={n_neighbors})",
            cfg.plot_dir / "variant_umap_by_n_ancestries_exact.png",
            legend_title="# ancestries", labels=labels)

    # Save a summary table for interpretation
    out = vdf.sort_values("mean_abs_weight", ascending=False).head(500).copy()
    out.to_csv(cfg.plot_dir / "top_variants_with_primary_ancestry_trait_sharedness.csv", index=False)

    print("Saved variant plots to:", cfg.plot_dir.resolve())
    print("Also saved:", (cfg.plot_dir / "top_variants_with_primary_ancestry_trait_sharedness.csv").name)


if __name__ == "__main__":
    main()
