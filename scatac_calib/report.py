"""统一输出：JSON + 最低报告清单 Markdown。"""

import json
from pathlib import Path

CHECKLIST = [
    ("zero_model", "零模型 AUROC（逐峰/逐细胞）", "zero"),
    ("delta_vs_zero", "ΔAUROC 与配对 bootstrap CI/p", "delta"),
    ("oracle", "群级 oracle（逐峰/逐细胞）", "oracle"),
    ("depth_matched", "深度分层后的逐峰 AUROC", "depth"),
    ("ap", "平均精度 AP", "ap"),
    ("m1", "M1 分半信度（有同源组时）", "m1"),
    ("mapq", "MAPQ≤1 比例 / 唯一比对比例（有 BAM 时）", "mapq"),
    ("random_peaks", "随机峰对照（可选）", "random"),
]


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)
    return path


def write_checklist(path, result, dataset_name=""):
    lines = [f"# 最低报告清单：{dataset_name or '未命名数据集'}", ""]
    lines.append("| 项目 | 说明 | 状态 |")
    lines.append("|---|---|---|")
    for key, desc, flag in CHECKLIST:
        present = _has(result, flag)
        lines.append(f"| {key} | {desc} | {'✅' if present else '— 缺失'} |")
    lines += ["", "## 主要数字", ""]
    for k, v in _flatten_numbers(result).items():
        lines.append(f"- **{k}**: {v}")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


def _has(result, flag):
    if flag == "zero":
        return "zero_model" in result
    if flag == "delta":
        return "delta_vs_zero" in result
    if flag == "oracle":
        return "oracle" in result
    if flag == "depth":
        return "depth_matched" in result
    if flag == "ap":
        return "ap" in result
    return flag in result


def _flatten_numbers(result, prefix=""):
    out = {}
    for k, v in result.items():
        name = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(_flatten_numbers(v, f"{name}."))
        elif isinstance(v, (int, float)):
            out[name] = v
    return out
