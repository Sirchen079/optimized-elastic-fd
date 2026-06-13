"""
Multi-panel bar chart: Taylor vs optimized coefficients quantitative error
against finer-grid reference, by snapshot time.

Reads:  results_standard/exp3_quant_error_vs_reference.csv
Writes: results_standard/figures/exp3_quant_error_comparison.png
        results_standard/figures/exp3_quant_error_comparison.pdf
"""

from pathlib import Path
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams


ROOT = Path(__file__).resolve().parent.parent / "results_standard"
CSV_FILE = ROOT / "exp3_quant_error_vs_reference.csv"
OUT_DIR = ROOT / "figures"
OUT_PNG = OUT_DIR / "exp3_quant_error_comparison.png"
OUT_PDF = OUT_DIR / "exp3_quant_error_comparison.pdf"


def configure_chinese_font():
    for nm in ["Microsoft YaHei", "SimHei", "SimSun", "FangSong", "KaiTi"]:
        rcParams["font.sans-serif"] = [nm] + rcParams["font.sans-serif"]
        rcParams["axes.unicode_minus"] = False
        return nm
    return None


def read_rows(path):
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fmt_value(v, scale_span):
    abs_v = abs(v)
    if abs_v >= 100:
        return f"{v:.2f}"
    if abs_v >= 1:
        return f"{v:.3f}"
    if abs_v >= 0.01 or scale_span >= 0.01:
        return f"{v:.4f}"
    return f"{v:.2e}"


def compute_ylim(vmin, vmax, span):
    if vmin >= 0:
        # all-positive: zoom in if range is narrow relative to magnitude
        narrow = vmin > 0 and span < vmin * 0.5
        if narrow:
            return (vmin - 0.20 * span, vmax + 0.45 * span)
        return (0.0, vmax + 0.40 * span)
    if vmax <= 0:
        # all-negative: bars hang from 0; leave headroom ABOVE 0 for badge
        return (vmin - 0.12 * span, 0 + 0.40 * span)
    return (vmin - 0.10 * span, vmax + 0.40 * span)


def main():
    configure_chinese_font()
    rows = read_rows(CSV_FILE)
    times = [float(r["time_s"]) for r in rows]
    n = len(rows)

    # name, csv-key-T, csv-key-O, panel-title, unit, higher_is_better, mode
    panels = [
        ("aligned_rel_L2", "taylor_aligned_rel_l2", "optimized_aligned_rel_l2",
         "对齐相对 L2 误差", "—", False, "pct"),
        ("RMSE", "taylor_rmse", "optimized_rmse",
         "均方根误差 RMSE", "m/s", False, "pct"),
        ("SNR(dB)", "taylor_snr_db", "optimized_snr_db",
         "信噪比 SNR", "dB", True, "delta_db"),
        ("Pearson", "taylor_pearson", "optimized_pearson",
         "Pearson 相关系数", "—", True, "delta"),
    ]

    color_taylor = "#3B6FB6"
    color_opt    = "#D7642C"
    color_good   = "#1a8917"
    color_bad    = "#c62828"
    bar_width = 0.36

    n_panels = len(panels)
    n_cols = 2
    n_rows = (n_panels + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12.5, 9.6))
    axes = axes.flatten()

    x = np.arange(n)
    tlabels = [f"t = {t:.3f} s" for t in times]

    for idx, (name, k_t, k_o, title, unit, higher_better, mode) in enumerate(panels):
        ax = axes[idx]
        vals_t = np.array([float(r[k_t]) for r in rows])
        vals_o = np.array([float(r[k_o]) for r in rows])

        b1 = ax.bar(x - bar_width / 2, vals_t, bar_width,
                    color=color_taylor, edgecolor="black", linewidth=0.6,
                    zorder=3)
        b2 = ax.bar(x + bar_width / 2, vals_o, bar_width,
                    color=color_opt, edgecolor="black", linewidth=0.6,
                    zorder=3)

        vmax = float(max(vals_t.max(), vals_o.max()))
        vmin = float(min(vals_t.min(), vals_o.min()))
        span = vmax - vmin if vmax > vmin else max(abs(vmax), 1e-12)

        ylo, yhi = compute_ylim(vmin, vmax, span)
        ax.set_ylim(ylo, yhi)
        plot_span = yhi - ylo

        all_negative = vmax <= 0

        # —— 柱内数值（白字，加描边以确保对比度）
        from matplotlib import patheffects as pe
        white_with_outline = [pe.withStroke(linewidth=1.6, foreground="black")]
        for bar in list(b1) + list(b2):
            h = bar.get_height()
            if h >= 0:
                y = h - 0.040 * plot_span
                va = "top"
            else:
                y = h + 0.040 * plot_span
                va = "bottom"
            ax.text(bar.get_x() + bar.get_width() / 2, y,
                    fmt_value(h, span),
                    ha="center", va=va, fontsize=9.5, color="white",
                    fontweight="bold", zorder=5,
                    path_effects=white_with_outline)

        # —— 柱顶上方“改善百分比 / Δ”徽标
        for i, (vt, vo) in enumerate(zip(vals_t, vals_o)):
            if mode == "pct":
                if vt != 0:
                    pct = (vt - vo) / vt * 100.0
                    tag = f"{pct:+.2f} %"
                else:
                    pct, tag = 0.0, "—"
                good = pct > 0
            elif mode == "delta_db":
                d = vo - vt
                tag = f"Δ = {d:+.3f} dB"
                good = d > 0
            else:  # delta
                d = vo - vt
                tag = f"Δ = {d:+.4f}"
                good = d > 0
            edge = color_good if good else color_bad

            if all_negative:
                # 把徽标固定在 y=0 与 yhi 之间的空白区（保证在轴范围内，避免被裁剪）
                badge_y = yhi - 0.55 * (yhi - 0.0)
            else:
                badge_y = max(vt, vo) + 0.18 * plot_span
                badge_y = min(badge_y, yhi - 0.06 * plot_span)

            ax.text(i, badge_y, tag,
                    ha="center", va="center",
                    fontsize=10, fontweight="bold", color=edge,
                    bbox=dict(boxstyle="round,pad=0.32",
                              facecolor="white", edgecolor=edge,
                              linewidth=1.0, alpha=0.95),
                    zorder=6)

        # —— 对负值条形图，画一条 y=0 的参考线，让“柱从 0 起”更清晰
        if all_negative or (vmin < 0 < vmax):
            ax.axhline(0, color="#888", linewidth=0.7, linestyle="-", zorder=1)

        ax.set_xticks(x)
        ax.set_xticklabels(tlabels, fontsize=10.5)
        ax.set_title(title, fontsize=12, pad=8)
        ax.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.6, zorder=0)
        ax.set_axisbelow(True)
        if unit and unit != "—":
            ax.set_ylabel(f"({unit})", fontsize=10)

    # 删除多余子图
    for j in range(n_panels, len(axes)):
        fig.delaxes(axes[j])

    # —— 标题 + 图例（两行布局，互不遮挡）
    fig.suptitle("exp3：优化系数 vs Taylor 系数 — 对细网格参考解的逐快照定量误差",
                 fontsize=14, y=0.985, fontweight="bold")

    handles = [plt.Rectangle((0, 0), 1, 1, color=color_taylor, ec="black"),
               plt.Rectangle((0, 0), 1, 1, color=color_opt, ec="black")]
    fig.legend(handles, ["Taylor 系数", "优化系数"],
               loc="upper center", ncol=2, fontsize=11.5,
               bbox_to_anchor=(0.5, 0.945), frameon=False,
               handlelength=1.6, handletextpad=0.6, columnspacing=2.0)

    fig.tight_layout(rect=[0, 0, 1, 0.92])
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
    fig.savefig(OUT_PDF, bbox_inches="tight")
    plt.close(fig)
    print("Wrote:")
    print("  ", OUT_PNG)
    print("  ", OUT_PDF)


if __name__ == "__main__":
    main()
