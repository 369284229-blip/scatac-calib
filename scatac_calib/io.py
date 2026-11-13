"""输入/输出：矩阵、预测、分组。

矩阵约定：所有内部计算都用 **峰 × 细胞**（n_peaks, n_cells）。输入可以显式
声明方向（``--matrix-orientation``）：``.npz`` 默认按峰 × 细胞读取；
``.h5ad`` 遵循 AnnData 的 ``obs × var`` 约定（细胞 × 峰），``auto`` 会先看
obs/var 的名字像不像基因组区间，再决定是否转置。
"""

import re
from pathlib import Path

import numpy as np
from scipy import sparse

# chr3:1000-2000 / chr3_1000_2000 / 3:1000-2000 / chr3|1000|2000
_PEAK_NAME_RE = re.compile(
    r"^(?:chr)?[A-Za-z0-9_.]+(?:[:_|-]\d+){2}$", re.IGNORECASE
)

ORIENTATIONS = ("auto", "peaks-x-cells", "cells-x-peaks")


def _fraction_peak_names(names):
    names = list(names)
    if not names:
        return 0.0
    hits = sum(1 for n in names if _PEAK_NAME_RE.match(str(n)))
    return hits / len(names)


def infer_h5ad_orientation(obs_names, var_names):
    """推断 AnnData 的方向，返回 ``peaks-x-cells`` 或 ``cells-x-peaks``。

    优先看名字：var 像基因组区间、obs 不像 -> 细胞 × 峰；反之 -> 峰 × 细胞。
    都像或都不像时退回 AnnData 的默认约定（观测 × 变量 = 细胞 × 峰）。
    """
    obs_like_peaks = _fraction_peak_names(obs_names)
    var_like_peaks = _fraction_peak_names(var_names)
    if var_like_peaks >= 0.5 and obs_like_peaks < 0.5:
        return "cells-x-peaks"
    if obs_like_peaks >= 0.5 and var_like_peaks < 0.5:
        return "peaks-x-cells"
    return "cells-x-peaks"


def _apply_orientation(matrix, source_orientation):
    """把 ``source_orientation`` 的矩阵转成峰 × 细胞。"""
    if source_orientation == "cells-x-peaks":
        return matrix.T.tocsr()
    return matrix


def load_matrix(path, orientation="auto", return_info=False):
    """加载矩阵并统一成 CSR (n_peaks, n_cells)。

    参数
    ----
    path : str | Path
        ``.npz``（scipy CSR 四件套）或 ``.h5ad``（需 scanpy）。
    orientation : {"auto", "peaks-x-cells", "cells-x-peaks"}
        ``.npz`` 无元数据可判，``auto`` 按峰 × 细胞读取；``.h5ad`` 的 ``auto``
        用 :func:`infer_h5ad_orientation` 判断。显式取值会校验并与推断结果对照。
    return_info : bool
        为 True 时返回 ``(matrix, info)``；``info`` 记录实际使用的方向。
    """
    if orientation not in ORIENTATIONS:
        raise ValueError(
            f"unknown matrix orientation {orientation!r}; expected one of {ORIENTATIONS}"
        )
    path = str(path)
    info = {"input": path, "format": None, "orientation": None,
            "transposed": False}
    if path.endswith(".npz"):
        d = np.load(path, allow_pickle=True)
        m = sparse.csr_matrix(
            (d["data"], d["indices"], d["indptr"]), shape=tuple(d["shape"])
        )
        m = sparse.csr_matrix(m)
        # .npz 没有名字可推断：auto 视为已按峰 × 细胞存储。
        stored = "peaks-x-cells"
        chosen = stored if orientation == "auto" else orientation
        info["format"] = "npz"
        info["inferred"] = False
        m = _apply_orientation(m, chosen)
        info["orientation"] = chosen
        info["transposed"] = chosen == "cells-x-peaks"
    elif path.endswith(".h5ad"):
        try:
            import scanpy as sc
        except ImportError as exc:  # pragma: no cover - exercised only without scanpy
            raise ImportError(
                "reading .h5ad requires scanpy; install with "
                "`pip install 'scatac-calib[h5ad]'`"
            ) from exc
        ad = sc.read_h5ad(path)
        inferred = infer_h5ad_orientation(ad.obs_names, ad.var_names)
        info["format"] = "h5ad"
        info["inferred"] = orientation == "auto"
        chosen = inferred if orientation == "auto" else orientation
        if orientation != "auto" and orientation != inferred:
            # 显式覆盖是允许的：用户知道文件不合惯例。
            info["inferred_orientation"] = inferred
        m = ad.X if sparse.issparse(ad.X) else sparse.csr_matrix(ad.X)
        m = sparse.csr_matrix(m)
        m = _apply_orientation(m, chosen)
        info["orientation"] = chosen
        info["transposed"] = chosen == "cells-x-peaks"
        info["n_obs"] = int(ad.n_obs)
        info["n_vars"] = int(ad.n_vars)
    else:
        raise ValueError(f"unsupported matrix format: {path}")
    m = m.astype(np.float32)
    info["shape"] = [int(m.shape[0]), int(m.shape[1])]
    if return_info:
        return m, info
    return m


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
