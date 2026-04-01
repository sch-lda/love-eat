#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PCoA workflow equivalent to the provided R script.

Features:
1. Read metadata and OTU abundance from Excel.
2. Align sample order by metadata Site column.
3. Compute Bray-Curtis distances across samples.
4. Run classical PCoA (cmdscale style).
5. Plot group points with adaptive shape (ellipse or line+triangle).
6. Save figure to PCoA_Result.png.
"""

from __future__ import annotations

import importlib
import subprocess
import sys
import warnings
from pathlib import Path


def ensure_packages() -> None:
    required = {
        "numpy": "numpy",
        "pandas": "pandas",
        "scipy": "scipy",
        "matplotlib": "matplotlib",
        "openpyxl": "openpyxl",
    }
    missing = []
    for pkg_name, import_name in required.items():
        try:
            importlib.import_module(import_name)
        except ImportError:
            missing.append(pkg_name)

    if missing:
        print(f"检测到缺失依赖，正在安装: {', '.join(missing)}")
        subprocess.check_call([sys.executable, "-m", "pip", "install", *missing])


ensure_packages()

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from scipy.spatial.distance import pdist, squareform


def pcoa_cmdscale(distance_matrix: np.ndarray, n_components: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """Classical PCoA (cmdscale) from a square distance matrix."""
    n = distance_matrix.shape[0]
    if distance_matrix.shape[0] != distance_matrix.shape[1]:
        raise ValueError("distance_matrix 必须是方阵。")

    d2 = distance_matrix ** 2
    i = np.eye(n)
    one = np.ones((n, n)) / n
    j = i - one
    b = -0.5 * j @ d2 @ j

    eigvals, eigvecs = np.linalg.eigh(b)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    pos = eigvals > 0
    eigvals_pos = eigvals[pos]
    eigvecs_pos = eigvecs[:, pos]

    if len(eigvals_pos) == 0:
        raise ValueError("PCoA 失败：没有正特征值。")

    coords = eigvecs_pos[:, :n_components] * np.sqrt(eigvals_pos[:n_components])
    return coords, eigvals


def add_confidence_ellipse(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    color: tuple[float, float, float, float],
    alpha: float = 0.12,
    n_std: float = 2.0,
    lw: float = 1.5,
) -> None:
    """Draw a covariance-based confidence ellipse for one group."""
    if x.size < 2 or y.size < 2:
        return

    cov = np.cov(x, y)
    if np.any(~np.isfinite(cov)):
        return

    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals = vals[order]
    vecs = vecs[:, order]

    if np.any(vals <= 0):
        return

    angle = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    width, height = 2 * n_std * np.sqrt(vals)
    center = (np.mean(x), np.mean(y))

    ell = Ellipse(
        xy=center,
        width=width,
        height=height,
        angle=angle,
        facecolor=color,
        edgecolor=color,
        lw=lw,
        alpha=alpha,
    )
    ax.add_patch(ell)


def is_group_near_linear(x: np.ndarray, y: np.ndarray, ratio_threshold: float = 0.08) -> bool:
    """Return True if group is too small or covariance is close to a line."""
    if x.size < 3 or y.size < 3:
        return True

    cov = np.cov(x, y)
    if np.any(~np.isfinite(cov)):
        return True

    vals = np.linalg.eigvalsh(cov)
    vals = np.sort(vals)
    if vals[-1] <= 0:
        return True

    # Minor-to-major axis variance ratio; smaller means more line-like.
    print(f"Group linearity ratio: {vals[0] / vals[-1]:.6f}")
    return (vals[0] / vals[-1]) < ratio_threshold


def order_points_along_major_axis(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Sort points along first principal direction for cleaner connecting lines."""
    pts = np.column_stack([x, y])
    centered = pts - pts.mean(axis=0)

    if pts.shape[0] < 2:
        return np.arange(pts.shape[0])

    cov = np.cov(centered, rowvar=False)
    vals, vecs = np.linalg.eigh(cov)
    major_vec = vecs[:, np.argmax(vals)]
    proj = centered @ major_vec
    return np.argsort(proj)


def main() -> None:
    meta_file = Path("./h1.xlsx")
    otu_file = Path("./h2.xlsx")
    out_file = Path("PCoA_Result.png")

    # PCoA 轴方向可按参考图调整。
    # 说明：PCoA 特征向量符号本身不唯一，不同软件实现常出现左右/上下镜像。
    flip_pcoa1 = False
    flip_pcoa2 = False

    print("正在读取数据...")
    meta_data = pd.read_excel(meta_file, sheet_name=0)
    otu_data = pd.read_excel(otu_file, sheet_name=0)

    required_meta_cols = {"Site", "Group"}
    if not required_meta_cols.issubset(meta_data.columns):
        raise ValueError(f"元数据文件缺少必要列: {required_meta_cols}")

    # 第一列视为 Taxonomy，其余列为样本丰度
    otu_matrix = otu_data.iloc[:, 1:].copy()
    otu_matrix.index = otu_data.iloc[:, 0].astype(str)

    sample_ids = meta_data["Site"].astype(str).tolist()
    otu_cols = otu_matrix.columns.astype(str).tolist()

    missing = sorted(set(sample_ids) - set(otu_cols))
    if missing:
        warnings.warn(f"以下样本在OTU表中未找到: {', '.join(missing)}")

    present_samples = [s for s in sample_ids if s in otu_cols]
    if len(present_samples) < 2:
        raise ValueError("可用于分析的样本少于2个，无法计算距离。")

    otu_sorted = otu_matrix[present_samples].apply(pd.to_numeric, errors="coerce").fillna(0.0)

    print("正在计算Bray-Curtis距离...")
    # 行=物种，列=样本；计算样本间距离需转置为 行=样本
    sample_matrix = otu_sorted.T.to_numpy(dtype=float)
    dist_vec = pdist(sample_matrix, metric="braycurtis")
    dist_mat = squareform(dist_vec)

    print("正在执行PCoA分析...")
    coords, eigvals_all = pcoa_cmdscale(dist_mat, n_components=2)

    pos_eig = eigvals_all[eigvals_all > 0]
    if len(pos_eig) < 2:
        raise ValueError("PCoA结果不足2个正特征值，无法绘制二维图。")

    pct = np.round(pos_eig[:2] / pos_eig.sum() * 100, 2)

    pcoa1 = -coords[:, 0] if flip_pcoa1 else coords[:, 0]
    pcoa2 = -coords[:, 1] if flip_pcoa2 else coords[:, 1]

    plot_df = pd.DataFrame(
        {
            "Sample": present_samples,
            "PCoA1": pcoa1,
            "PCoA2": pcoa2,
        }
    )

    final_df = plot_df.merge(
        meta_data[["Site", "Group"]].rename(columns={"Site": "Sample"}),
        on="Sample",
        how="left",
    )

    print("正在绘制PCoA图...")
    # 稍微调大画布比例，避免标签和图例拥挤
    fig, ax = plt.subplots(figsize=(11, 8), dpi=130)

    groups = final_df["Group"].astype(str).fillna("NA")
    unique_groups = groups.unique().tolist()
    cmap = plt.get_cmap("Set3", max(3, len(unique_groups)))
    color_map = {g: cmap(i) for i, g in enumerate(unique_groups)}

    for g in unique_groups:
        sub = final_df[groups == g]
        x = sub["PCoA1"].to_numpy()
        y = sub["PCoA2"].to_numpy()
        c = color_map[g]

        if is_group_near_linear(x, y, ratio_threshold=0.003754):
            order = order_points_along_major_axis(x, y)
            x_line = x[order]
            y_line = y[order]

            ax.plot(x_line, y_line, color=c, linewidth=1.4, alpha=0.9)
            ax.scatter(
                x,
                y,
                marker="^",
                s=62,
                alpha=0.92,
                color=c,
                edgecolors="white",
                linewidths=0.6,
                label=g,
            )
        else:
            add_confidence_ellipse(ax, x, y, color=c, alpha=0.14, n_std=2.0, lw=1.2)
            ax.scatter(x, y, s=45, alpha=0.9, color=c, edgecolors="white", linewidths=0.5, label=g)

    x_min, x_max = final_df["PCoA1"].min(), final_df["PCoA1"].max()
    y_min, y_max = final_df["PCoA2"].min(), final_df["PCoA2"].max()
    pad_x = max((x_max - x_min) * 0.08, 0.02)
    pad_y = max((y_max - y_min) * 0.08, 0.02)

    ax.set_xlim(x_min - pad_x, x_max + pad_x)
    ax.set_ylim(y_min - pad_y, y_max + pad_y)

    ax.axvline(0, linestyle="--", color="gray", alpha=0.6, linewidth=1.0)
    ax.axhline(0, linestyle="--", color="gray", alpha=0.6, linewidth=1.0)

    ax.set_xlabel(f"PCoA 1 ({pct[0]}%)", fontsize=12)
    ax.set_ylabel(f"PCoA 2 ({pct[1]}%)", fontsize=12)
    ax.set_title("PCoA Plot (Bray-Curtis Distance)", fontsize=14)

    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
        spine.set_color("black")

    ax.tick_params(axis="both", labelsize=10)
    ax.legend(title="Group", loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)

    fig.tight_layout()
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"分析完成！图片已保存为 {out_file}")


if __name__ == "__main__":
    main()
