# 03_train_eval_embed_svd.py
from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Tuple, Dict, List

import numpy as np
import pandas as pd
import scipy.sparse as sp

import torch
import torch.nn as nn
import torch.optim as optim

from sklearn.decomposition import TruncatedSVD


@dataclass
class Config:
    # Inputs
    x_csr_npz: Path = Path("./data/tensors/X_pgs_feat_csr.npz")
    pgs_ids_json: Path = Path("./data/tensors/pgs_ids.json")
    locus_index_json: Path = Path("./data/tensors/locus_index.json")
    feature_schema_json: Path = Path("./data/tensors/feature_schema.json")
    meta_parquet: Path = Path("./data/processed/pgs_metadata.parquet")

    # Outputs
    embedding_dir: Path = Path("./results/embeddings")

    # Split
    test_frac: float = 0.2
    seed: int = 42

    # Sparse compression (key change)
    svd_dim: int = 1024  # try 512/1024/2048 depending on #PGS
    svd_random_state: int = 42

    # AE on compressed inputs
    pgs_latent_dim: int = 32
    pgs_hidden_dim: int = 256
    pgs_epochs: int = 300
    lr: float = 1e-3
    weight_decay: float = 1e-5

    # Variant embeddings (weights-only channel)
    # With small #PGS, keep latent small
    var_latent_dim: int = 4
    var_hidden_dim: int = 64
    var_epochs: int = 200
    min_support: int = 2


def ensure_dirs(cfg: Config):
    cfg.embedding_dir.mkdir(parents=True, exist_ok=True)


def set_seed(seed: int):
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def device() -> torch.device:
    return torch.device("mps" if torch.backends.mps.is_available() else ("cuda" if torch.cuda.is_available() else "cpu"))


def train_test_split_indices(n: int, test_frac: float, seed: int) -> Tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    idx = np.arange(n)
    rng.shuffle(idx)
    n_test = max(1, int(round(n * test_frac)))
    return idx[n_test:], idx[:n_test]


class AE(nn.Module):
    def __init__(self, input_dim: int, hidden_dim: int, latent_dim: int, dropout: float = 0.1):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, latent_dim),
        )
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, input_dim),
        )

    def forward(self, x: torch.Tensor):
        z = self.encoder(x)
        x_hat = self.decoder(z)
        return x_hat, z


@torch.no_grad()
def reconstruction_metrics(x: torch.Tensor, x_hat: torch.Tensor) -> Dict[str, float]:
    mse = torch.mean((x_hat - x) ** 2).item()
    eps = 1e-8
    num = torch.sum(x * x_hat, dim=1)
    den = torch.norm(x, dim=1) * torch.norm(x_hat, dim=1) + eps
    cos = torch.mean(num / den).item()
    return {"mse": mse, "cosine": cos}


def main():
    cfg = Config()
    ensure_dirs(cfg)
    set_seed(cfg.seed)

    dev = device()
    print("Using device:", dev)
    print("Config:", json.dumps({k: str(v) for k, v in asdict(cfg).items()}, indent=2))

    X_csr = sp.load_npz(cfg.x_csr_npz).tocsr()
    pgs_ids: List[str] = json.loads(cfg.pgs_ids_json.read_text())
    locus_index = json.loads(cfg.locus_index_json.read_text())
    feat_schema = json.loads(cfg.feature_schema_json.read_text())
    meta = pd.read_parquet(cfg.meta_parquet)

    B, D = X_csr.shape
    print(f"Loaded CSR: B={B} PGS, D={D} features, nnz={X_csr.nnz}")
    assert B == len(pgs_ids)

    # ---------------- Step 1: Sparse compression (NO densify of X) ----------------
    # Note: svd_dim must be < B in many useful regimes, but sklearn allows it.
    svd_dim = min(cfg.svd_dim, max(2, B - 1)) if B <= cfg.svd_dim else cfg.svd_dim
    print(f"Running TruncatedSVD to {svd_dim} dims...")
    svd = TruncatedSVD(n_components=svd_dim, random_state=cfg.svd_random_state)
    X_svd = svd.fit_transform(X_csr).astype(np.float32)  # dense [B, svd_dim]
    explained = float(np.sum(svd.explained_variance_ratio_))
    (cfg.embedding_dir / "svd_explained_variance.txt").write_text(f"{explained}\n")
    print(f"SVD done. Explained variance ratio sum: {explained:.4f}")

    # L2 normalize rows (good for cosine geometry)
    X_svd = X_svd / (np.linalg.norm(X_svd, axis=1, keepdims=True) + 1e-8)

    # Save compressed features for plotting/debugging
    np.save(cfg.embedding_dir / "X_svd.npy", X_svd)
    torch.save({"pgs_ids": pgs_ids, "X_svd": torch.tensor(X_svd)}, cfg.embedding_dir / "X_svd.pt")

    # ---------------- Step 2: Train/test split on compressed X ----------------
    train_idx, test_idx = train_test_split_indices(B, cfg.test_frac, cfg.seed)
    X_train = torch.tensor(X_svd[train_idx], dtype=torch.float32).to(dev)
    X_test = torch.tensor(X_svd[test_idx], dtype=torch.float32).to(dev)
    X_all = torch.tensor(X_svd, dtype=torch.float32).to(dev)

    # ---------------- Step 3: AE on compressed features ----------------
    model = AE(input_dim=svd_dim, hidden_dim=cfg.pgs_hidden_dim, latent_dim=cfg.pgs_latent_dim).to(dev)
    opt = optim.Adam(model.parameters(), lr=cfg.lr, weight_decay=cfg.weight_decay)
    mse = nn.MSELoss()

    history = []
    for epoch in range(1, cfg.pgs_epochs + 1):
        model.train()
        opt.zero_grad()
        x_hat, _ = model(X_train)
        loss = mse(x_hat, X_train)
        loss.backward()
        opt.step()

        if epoch == 1 or epoch % 25 == 0:
            model.eval()
            with torch.no_grad():
                train_hat, _ = model(X_train)
                test_hat, _ = model(X_test)
            mt_train = reconstruction_metrics(X_train, train_hat)
            mt_test = reconstruction_metrics(X_test, test_hat)
            history.append({"epoch": epoch, "train": mt_train, "test": mt_test})
            print(f"[PGS-AE] epoch {epoch:03d} loss={loss.item():.6f} "
                  f"train_mse={mt_train['mse']:.6f} test_mse={mt_test['mse']:.6f}")

    model.eval()
    with torch.no_grad():
        _, Z_all = model(X_all)
    Z_all = Z_all.detach().cpu()

    torch.save(
        {
            "pgs_ids": pgs_ids,
            "embeddings": Z_all,
            "config": {
                "svd_dim": svd_dim,
                "svd_explained_variance_sum": explained,
                "pgs_latent_dim": cfg.pgs_latent_dim,
                "pgs_hidden_dim": cfg.pgs_hidden_dim,
                "pgs_epochs": cfg.pgs_epochs,
            },
        },
        cfg.embedding_dir / "pgs_embeddings.pt",
    )
    (cfg.embedding_dir / "pgs_ae_metrics.json").write_text(json.dumps({"history": history}, indent=2))

    # ---------------- Step 4: Variant embeddings (weights-only channel) ----------------
    # This part does not densify the full multi-feature matrix, only the weight_scaled slice.
    F = int(feat_schema["F"])
    wfi = int(feat_schema["feat_index"]["weight_scaled"])
    V = len(locus_index)

    weight_cols = np.array([i * F + wfi for i in range(V)], dtype=np.int64)
    Xw = X_csr[:, weight_cols].tocsr()  # [B, V]
    Xw_dense = Xw.toarray().astype(np.float32)  # [B, V] -- this is manageable if V is not huge

    X_var_full = torch.tensor(Xw_dense.T, dtype=torch.float32)  # [V, B]
    support = torch.count_nonzero(X_var_full, dim=1).numpy()
    keep_idx = np.where(support >= cfg.min_support)[0]
    X_var = X_var_full[keep_idx].to(dev)

    print(f"[VAR] loci={V}, kept={len(keep_idx)} with min_support={cfg.min_support}")

    var_model = AE(input_dim=B, hidden_dim=cfg.var_hidden_dim, latent_dim=cfg.var_latent_dim).to(dev)
    var_opt = optim.Adam(var_model.parameters(), lr=cfg.lr, weight_decay=cfg.weight_decay)

    for epoch in range(1, cfg.var_epochs + 1):
        var_model.train()
        var_opt.zero_grad()
        xh, _ = var_model(X_var)
        loss_v = mse(xh, X_var)
        loss_v.backward()
        var_opt.step()
        if epoch == 1 or epoch % 50 == 0:
            print(f"[VAR-AE] epoch {epoch:03d} loss={loss_v.item():.6f}")

    var_model.eval()
    with torch.no_grad():
        _, Z_var = var_model(X_var)
    Z_var = Z_var.detach().cpu().numpy()

    # invert locus index
    inv_locus = [""] * V
    for lid, j in locus_index.items():
        inv_locus[int(j)] = lid

    nnz_across = np.count_nonzero(Xw_dense, axis=0)
    mean_abs = np.mean(np.abs(Xw_dense), axis=0)

    var_df = pd.DataFrame({
        "locus_id": [inv_locus[i] for i in keep_idx],
        "locus_index": keep_idx,
        "nnz_across_pgs": nnz_across[keep_idx],
        "mean_abs_weight": mean_abs[keep_idx],
    })
    for k in range(Z_var.shape[1]):
        var_df[f"z{k}"] = Z_var[:, k]
    var_df.to_csv(cfg.embedding_dir / "variant_embeddings.csv", index=False)

    snapshot = {
        "B": B,
        "D": D,
        "svd_dim": svd_dim,
        "svd_explained_variance_sum": explained,
        "F": F,
        "V": V,
        "variant_keep": int(len(keep_idx)),
    }
    (cfg.embedding_dir / "run_snapshot.json").write_text(json.dumps(snapshot, indent=2))

    print("Done. Saved to:", cfg.embedding_dir.resolve())
    print(" - pgs_embeddings.pt")
    print(" - pgs_ae_metrics.json")
    print(" - X_svd.npy (for plotting/debug)")
    print(" - variant_embeddings.csv")


if __name__ == "__main__":
    main()
