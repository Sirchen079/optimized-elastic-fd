"""
Quantitative error comparison of Taylor vs. optimized (minimax) staggered-grid
schemes against the finer-grid reference solution from exp3.

For each snapshot, compute on the downsampled vz field:
  - raw relative L2 error                  ‖x - r‖₂ / ‖r‖₂
  - scale-aligned relative L2 error        ‖x - α r‖₂ / ‖α r‖₂   (α = <x,r>/<r,r>)
  - max absolute error                      max|x - r|
  - RMSE                                    sqrt(mean((x - r)^2))
  - NRMSE (vs reference range)              RMSE / (max(r) - min(r))
  - SNR in dB                               20·log10(‖r‖₂ / ‖x - r‖₂)
  - Pearson correlation coefficient         corr(x, r)

Outputs:
  results_standard/exp3_quant_error_vs_reference.csv
  results_standard/exp3_quant_error_summary.txt
"""

import sys
from pathlib import Path

import numpy as np
import scipy.io as sio


ROOT = Path(__file__).resolve().parent.parent / "results_standard"
REF_FILE = ROOT / "exp3_reference_finer_grid.mat"
TAYLOR_FILE = ROOT / "exp3_staggered_pml_standard.mat"
OPT_FILE = ROOT / "exp3_staggered_pml_minimax.mat"

OUT_CSV = ROOT / "exp3_quant_error_vs_reference.csv"
OUT_TXT = ROOT / "exp3_quant_error_summary.txt"


def load_result(path: Path):
    d = sio.loadmat(str(path))
    r = d["result"][0, 0]
    snaps = r["snapshots"][0, 0]
    vz = np.asarray(snaps["vz"], dtype=np.float64)
    t = np.asarray(snaps["time"], dtype=np.float64).flatten()
    dt = float(np.asarray(r["dt"]).flatten()[0])
    coeff = np.asarray(r["coeff"], dtype=np.float64).flatten()
    return {"vz": vz, "time": t, "dt": dt, "coeff": coeff,
            "tag": str(r["tag"][0])}


def downsample_to(coarse_shape, ref_field):
    nz_c, nx_c = coarse_shape
    nz_r, nx_r = ref_field.shape
    fz = round(nz_r / nz_c)
    fx = round(nx_r / nx_c)
    if fz < 1 or fx < 1:
        raise ValueError("Reference grid must be finer than candidate grid.")
    sub = ref_field[::fz, ::fx]
    return sub[:nz_c, :nx_c]


def metrics(candidate, reference):
    diff = candidate - reference
    norm_r = np.linalg.norm(reference)
    norm_c = np.linalg.norm(candidate)
    norm_d = np.linalg.norm(diff)

    # scale-aligned error (optimal scalar α)
    denom = float(np.dot(reference.flatten(), reference.flatten())) + np.finfo(float).eps
    alpha = float(np.dot(candidate.flatten(), reference.flatten())) / denom
    aligned_diff = candidate - alpha * reference
    aligned_err = float(np.linalg.norm(aligned_diff)) / (float(np.linalg.norm(alpha * reference)) + np.finfo(float).eps)

    rmse = float(np.sqrt(np.mean(diff ** 2)))
    rng = float(reference.max() - reference.min())
    nrmse = rmse / rng if rng > 0 else float("nan")

    rel_l2 = float(norm_d) / (float(norm_r) + np.finfo(float).eps)
    snr_db = 20.0 * np.log10((float(norm_r) + np.finfo(float).eps) /
                             (float(norm_d) + np.finfo(float).eps))
    max_abs = float(np.max(np.abs(diff)))

    # Pearson correlation
    c0 = candidate - candidate.mean()
    r0 = reference - reference.mean()
    denom_p = float(np.linalg.norm(c0) * np.linalg.norm(r0)) + np.finfo(float).eps
    pearson = float(np.dot(c0.flatten(), r0.flatten()) / denom_p)

    return {
        "rel_l2": rel_l2,
        "aligned_rel_l2": aligned_err,
        "alpha": alpha,
        "max_abs": max_abs,
        "rmse": rmse,
        "nrmse": nrmse,
        "snr_db": float(snr_db),
        "pearson": pearson,
        "norm_ref": float(norm_r),
        "norm_cand": float(norm_c),
        "ref_peak_abs": float(np.max(np.abs(reference))),
    }


def main():
    ref = load_result(REF_FILE)
    tay = load_result(TAYLOR_FILE)
    opt = load_result(OPT_FILE)

    nSnap = min(tay["vz"].shape[2], opt["vz"].shape[2], ref["vz"].shape[2])
    coarse_shape = tay["vz"].shape[:2]

    rows = []
    for k in range(nSnap):
        ref_k = downsample_to(coarse_shape, ref["vz"][:, :, k])
        tay_k = tay["vz"][:, :, k]
        opt_k = opt["vz"][:, :, k]

        m_t = metrics(tay_k, ref_k)
        m_o = metrics(opt_k, ref_k)

        rows.append({
            "k": k + 1,
            "time": float(tay["time"][k]),
            "ref_time": float(ref["time"][k]),
            "ref_peak_abs": m_t["ref_peak_abs"],
            "T": m_t,
            "O": m_o,
        })

    # ── Write CSV ──────────────────────────────────────────────────────────────
    header = [
        "snapshot_index", "time_s", "reference_time_s", "reference_peak_abs",
        "taylor_rel_l2", "optimized_rel_l2", "rel_l2_reduction_pct",
        "taylor_aligned_rel_l2", "optimized_aligned_rel_l2", "aligned_reduction_pct",
        "taylor_max_abs", "optimized_max_abs", "max_abs_reduction_pct",
        "taylor_rmse", "optimized_rmse", "rmse_reduction_pct",
        "taylor_nrmse", "optimized_nrmse",
        "taylor_snr_db", "optimized_snr_db", "snr_gain_db",
        "taylor_pearson", "optimized_pearson",
        "taylor_alpha", "optimized_alpha",
    ]
    lines_csv = [",".join(header)]
    for r in rows:
        T, O = r["T"], r["O"]

        def pct(t, o):
            return (t - o) / t * 100.0 if t != 0 else float("nan")

        line = [
            str(r["k"]),
            f"{r['time']:.15g}",
            f"{r['ref_time']:.15g}",
            f"{r['ref_peak_abs']:.15g}",
            f"{T['rel_l2']:.15g}",
            f"{O['rel_l2']:.15g}",
            f"{pct(T['rel_l2'], O['rel_l2']):.6f}",
            f"{T['aligned_rel_l2']:.15g}",
            f"{O['aligned_rel_l2']:.15g}",
            f"{pct(T['aligned_rel_l2'], O['aligned_rel_l2']):.6f}",
            f"{T['max_abs']:.15g}",
            f"{O['max_abs']:.15g}",
            f"{pct(T['max_abs'], O['max_abs']):.6f}",
            f"{T['rmse']:.15g}",
            f"{O['rmse']:.15g}",
            f"{pct(T['rmse'], O['rmse']):.6f}",
            f"{T['nrmse']:.15g}",
            f"{O['nrmse']:.15g}",
            f"{T['snr_db']:.6f}",
            f"{O['snr_db']:.6f}",
            f"{O['snr_db'] - T['snr_db']:.6f}",
            f"{T['pearson']:.15g}",
            f"{O['pearson']:.15g}",
            f"{T['alpha']:.15g}",
            f"{O['alpha']:.15g}",
        ]
        lines_csv.append(",".join(line))

    OUT_CSV.write_text("\n".join(lines_csv) + "\n", encoding="utf-8")

    # ── Aggregate / mean values ────────────────────────────────────────────────
    def mean_of(key, side):
        return float(np.mean([r[side][key] for r in rows]))

    mean_T = {k: mean_of(k, "T") for k in
              ["rel_l2", "aligned_rel_l2", "max_abs", "rmse", "nrmse", "snr_db", "pearson"]}
    mean_O = {k: mean_of(k, "O") for k in
              ["rel_l2", "aligned_rel_l2", "max_abs", "rmse", "nrmse", "snr_db", "pearson"]}

    def pct_or_nan(t, o):
        return (t - o) / t * 100.0 if t != 0 else float("nan")

    # ── Write human-readable TXT ───────────────────────────────────────────────
    txt = []
    txt.append("===== exp3 优化系数 vs Taylor 系数 — 对细网格参考解的定量误差对比 =====")
    txt.append("")
    txt.append("数据来源：")
    txt.append(f"  参考解（细网格）       : {REF_FILE.name}")
    txt.append(f"     dt = {ref['dt']:.6e} s,  vz 形状 = {ref['vz'].shape}")
    txt.append(f"  Taylor（粗网格）       : {TAYLOR_FILE.name}")
    txt.append(f"     dt = {tay['dt']:.6e} s,  vz 形状 = {tay['vz'].shape}")
    txt.append(f"  优化系数（粗网格）     : {OPT_FILE.name}")
    txt.append(f"     dt = {opt['dt']:.6e} s,  vz 形状 = {opt['vz'].shape}")
    txt.append("")
    txt.append("空间对齐：参考解每隔 (round(nz_r/nz_c), round(nx_r/nx_c)) 抽样到粗网格。")
    txt.append("时间对齐：使用各方案保存的 snapshot 实际时刻（误差 < 4e-4 s，见 exp3_snapshot_time_alignment.csv）。")
    txt.append("说明：所有指标均在 vz 分量上计算（exp3 仅保存 vz 快照）。")
    txt.append("")
    txt.append("─── 指标定义 ──────────────────────────────────────────────────────")
    txt.append("  rel_L2          = ‖cand − ref‖₂ / ‖ref‖₂")
    txt.append("  aligned_rel_L2  = ‖cand − α·ref‖₂ / ‖α·ref‖₂,  α = <cand,ref>/<ref,ref>")
    txt.append("  max_abs         = max|cand − ref|")
    txt.append("  RMSE            = sqrt( mean( (cand − ref)² ) )")
    txt.append("  NRMSE           = RMSE / (max(ref) − min(ref))")
    txt.append("  SNR(dB)         = 20·log10( ‖ref‖₂ / ‖cand − ref‖₂ )")
    txt.append("  Pearson         = ⟨cand−mean, ref−mean⟩ / (‖·‖·‖·‖)")
    txt.append("")

    for r in rows:
        T, O = r["T"], r["O"]
        txt.append(f"─── 快照 #{r['k']}  t ≈ {r['time']:.6f} s  (参考解 t = {r['ref_time']:.6f} s) ───")
        txt.append(f"  参考波场峰值幅值 = {r['ref_peak_abs']:.6e}")
        txt.append("")
        txt.append("  指标               Taylor               优化系数             相对降低")
        def line_of(name, key, fmt="{:>16.6e}"):
            t = T[key]; o = O[key]
            red = pct_or_nan(t, o)
            return f"  {name:<14} " + fmt.format(t) + "    " + fmt.format(o) + f"    {red:+8.3f} %"
        txt.append(line_of("rel_L2",         "rel_l2"))
        txt.append(line_of("aligned_rel_L2", "aligned_rel_l2"))
        txt.append(line_of("max_abs",        "max_abs"))
        txt.append(line_of("RMSE",           "rmse"))
        txt.append(line_of("NRMSE",          "nrmse"))
        # SNR / Pearson 没有"降低 %"概念,改为"增益 dB"
        snr_t = T["snr_db"]; snr_o = O["snr_db"]
        txt.append(f"  SNR(dB)        {snr_t:>16.6f}    {snr_o:>16.6f}    {snr_o - snr_t:+8.3f} dB (增益)")
        p_t = T["pearson"]; p_o = O["pearson"]
        txt.append(f"  Pearson        {p_t:>16.9f}    {p_o:>16.9f}    {p_o - p_t:+8.6f} (Δ)")
        txt.append(f"  α (能量缩放)    {T['alpha']:>16.9f}    {O['alpha']:>16.9f}")
        txt.append("")

    txt.append("─── 三个快照的均值 ─────────────────────────────────────────────────")
    txt.append("  指标               Taylor               优化系数             相对降低")
    for name, key in [
        ("rel_L2",         "rel_l2"),
        ("aligned_rel_L2", "aligned_rel_l2"),
        ("max_abs",        "max_abs"),
        ("RMSE",           "rmse"),
        ("NRMSE",          "nrmse"),
    ]:
        t = mean_T[key]; o = mean_O[key]
        red = pct_or_nan(t, o)
        txt.append(f"  {name:<14} {t:>16.6e}    {o:>16.6e}    {red:+8.3f} %")
    txt.append(f"  SNR(dB)        {mean_T['snr_db']:>16.6f}    {mean_O['snr_db']:>16.6f}    {mean_O['snr_db'] - mean_T['snr_db']:+8.3f} dB (增益)")
    txt.append(f"  Pearson        {mean_T['pearson']:>16.9f}    {mean_O['pearson']:>16.9f}    {mean_O['pearson'] - mean_T['pearson']:+8.6f} (Δ)")
    txt.append("")

    txt.append("─── 与既有 exp3_reference_error_metrics.txt 中数值的对应 ──────────")
    txt.append("  本文件中 aligned_rel_L2 ≡ exp3_reference_error_metrics.txt 中的")
    txt.append("  '相对参考波场偏离' 列（最优能量缩放后的归一化 L2 残差）。")
    txt.append("")

    OUT_TXT.write_text("\n".join(txt) + "\n", encoding="utf-8")

    # ── Echo summary to stdout ────────────────────────────────────────────────
    print("Wrote:")
    print("  ", OUT_CSV)
    print("  ", OUT_TXT)
    print()
    print("Mean over 3 snapshots:")
    for name, key in [
        ("rel_L2", "rel_l2"),
        ("aligned_rel_L2", "aligned_rel_l2"),
        ("max_abs", "max_abs"),
        ("RMSE", "rmse"),
        ("NRMSE", "nrmse"),
    ]:
        t = mean_T[key]; o = mean_O[key]
        red = pct_or_nan(t, o)
        print(f"  {name:<14} Taylor={t:.6e}  Optim={o:.6e}  reduction={red:+.3f}%")
    print(f"  SNR(dB)        Taylor={mean_T['snr_db']:.4f}  Optim={mean_O['snr_db']:.4f}  gain={mean_O['snr_db']-mean_T['snr_db']:+.4f} dB")
    print(f"  Pearson        Taylor={mean_T['pearson']:.9f}  Optim={mean_O['pearson']:.9f}")


if __name__ == "__main__":
    sys.exit(main())
