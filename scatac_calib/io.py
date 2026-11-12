"""输入/输出：矩阵、预测、分组。"""

from pathlib import Path

import numpy as np
from scipy import sparse


def load_matrix(path):
    """加载 cell x peak 矩阵。

    支持:
      - .npz: scipy CSR 四件套 (data/indices/indptr/shape)
      - .h5ad: AnnData（需 scanpy）
    返回 CSR (n_peaks, n_cells)。
    """
    path = str(path)
    if path.endswith(".npz"):
        d = np.load(path, allow_pickle=True)
        m = sparse.csr_matrix(
            (d["data"], d["indices"], d["indptr"]), shape=tuple(d["shape"])
        )
    elif path.endswith(".h5ad"):
        import scanpy as sc

        ad = sc.read_h5ad(path)
        m = ad.X if sparse.issparse(ad.X) else sparse.csr_matrix(ad.X)
        m = sparse.csr_matrix(m)
    else:
        raise ValueError(f"unsupported matrix format: {path}")
    return m.astype(np.float32)


def load_preds(path, key="preds"):
    """加载模型预测（dense .npy / .npz），形状 (n_peaks, n_cells)。"""
    path = str(path)
    if path.endswith(".npz"):
        d = np.load(path, allow_pickle=True)
        if key not in d:
            raise KeyError(f"{key!r} not in {path}; keys={list(d.keys())}")
        arr = d[key]
    else:
        arr = np.load(path)
    if sparse.issparse(arr):
        arr = arr.toarray()
    return np.asarray(arr, dtype=np.float32)


def load_labels(path, key="labels"):
    """从 preds.npz 这类文件里取标签矩阵。"""
    d = np.load(str(path), allow_pickle=True)
    if key not in d:
        raise KeyError(f"{key!r} not found; keys={list(d.keys())}")
    arr = d[key]
    if sparse.issparse(arr):
        arr = arr.toarray()
    return np.asarray(arr, dtype=np.float32)


def load_groups(path):
    """读取 cell_groups.txt（barcode<TAB>cluster）或一行一个标签。"""
    labels = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            labels.append(parts[1] if len(parts) > 1 else parts[0].split()[-1])
    return np.array(labels)


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)
    return Path(path)
