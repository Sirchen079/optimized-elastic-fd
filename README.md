# optimized-elastic-fd

> 优化差分系数的交错网格有限差分弹性波模拟 · Optimized Staggered-Grid Finite-Difference Elastic Wave Modeling

MATLAB 项目。基于 Marmousi 弹性模型，实现交错网格有限差分弹性波正演模拟（含 PML 吸收边界），并通过模拟退火算法以最大范数频散误差为目标优化差分系数，降低数值频散。

## 目录结构

```
.
├── run_ewm.m                  # 主入口脚本
├── src/                           # 全部 MATLAB 源代码
│   ├── ewm_staggered_solver.m     # 交错网格有限差分求解器（含 PML）
│   ├── ewm_regular_solver.m       # 常规网格求解器（仅实验1对比用）
│   ├── ewm_pml2d.m                # PML 完美匹配层实现
│   ├── ewm_optimize_minimax_coeffs.m  # 模拟退火优化差分系数
│   ├── ewm_fd_coefficients.m      # Taylor 展开差分系数计算
│   ├── ewm_dispersion_metrics.m   # 频散误差分析
│   ├── ewm_default_config.m       # 实验配置（preview/standard）
│   ├── ewm_load_marmousi.m        # Marmousi 模型加载与重采样
│   ├── ewm_prepare_elastic_params.m   # 弹性参数 λ, μ, ρ 预处理
│   ├── ewm_ricker.m               # Ricker 子波震源函数
│   ├── ewm_stable_dt.m            # CFL 稳定时间步长计算
│   ├── ewm_pad2d.m                # 二维边界填充
│   ├── ewm_wavefield_colormap.m   # 蓝-白-红发散色图（波场可视化）
│   ├── ewm_apply_chinese_style.m  # 中文字体设置
│   ├── ewm_save_figure.m          # 图形保存（300 DPI）
│   ├── ewm_plot_*.m               # 各类绘图函数（共 20+ 个）
│   ├── ewm_write_summary.m        # 结果摘要输出
│   ├── ewm_write_required_tables.m    # 参数表格输出
│   └── ...
├── scripts/                       # 辅助脚本（独立绘图、高 DPI 重绘、定量分析）
├── models/                        # 速度模型目录（放入自己的 Marmousi 模型，详见 models/README.md）
├── results_preview/                 # 快速验证模式输出（80m 网格）
│   └── figures/                   # 所有图形输出
├── results_standard/                # 论文模式输出（10m 网格）
│   └── figures/
└── results_const_velocity/        # 常速度模型实验输出
```

## 运行方式

```matlab
run_ewm            % 快速验证（80m 网格，约 2 分钟）
run_ewm('preview')   % 同上
run_ewm('standard')  % 论文模式（10m 网格，约 30 分钟）
```

两种模式的核心区别：

| 参数 | preview | standard |
|------|---------|----------|
| 网格间距 | 80 m | 10 m |
| Ricker 主频 | 4 Hz | 12 Hz |
| PML 层数 | 12 | 30 |
| 时间步数 | ~360 | ~1600 |

## 三个实验

### 实验 1：常规网格 vs 交错网格

对比常规网格与交错网格两种有限差分离散方案的波场差异。交错网格将速度分量和应力分量错开半个网格存放，空间精度更高，是弹性波模拟的主流方法。

- 快照时刻：1.00 s, 1.10 s, 1.20 s
- 无吸收边界，标准 Taylor 系数

### 实验 2：无吸收边界 vs PML 吸收边界

对比无吸收边界条件（人工截断）与 PML 完美匹配层吸收边界的效果。无吸收边界会产生明显的边界反射波，PML 则将出射波在边界层内逐步衰减。

- 快照时刻：1.65 s（此时边界反射已充分发展）
- 交错网格，标准 Taylor 系数

### 实验 3：标准 Taylor 系数 vs 基于最大范数目标函数的优化系数

对比传统 Taylor 展开差分系数与通过模拟退火优化的 minimax 差分系数。优化目标为最小化给定波数范围内的最大绝对频散误差 max|k_num·Δ - k·Δ|。

- 快照时刻：0.65 s, 0.95 s, 1.25 s
- 交错网格 + PML 吸收边界

---

## 输出图形详细说明

每次运行会在对应 `results_*/figures/` 目录下生成以下图形。

### 基础模型与震源

| 文件名 | 含义 |
|--------|------|
| `marmousi_model.png` | **Marmousi 弹性模型三参数面板图**。从左到右分别显示 P 波速度 Vp (km/s)、S 波速度 Vs (km/s)、密度 ρ (g/cm³) 的空间分布。横轴为水平距离，纵轴为深度（向下增加）。该模型是国际通用的复杂地质构造基准模型，包含褶皱、断层、速度反转等特征。|
| `vs_model.png` | **S 波速度单独大图**。更清晰地展示 Vs 分布细节及人工偏移量（若有）。|
| `ricker_wavelet_time_spectrum.png` | **Ricker 子波时域波形与频谱**。左图为归一化时域波形，右图为振幅谱。Ricker 子波是地震模拟中最常用的零相位震源信号，其频谱以主频 f₀ 为中心，带宽约 0 ~ 2.5f₀。|

### 参数与稳定性表格

| 文件名 | 含义 |
|--------|------|
| `experiment_parameters_table.png` | **完整实验参数表**。列出模型尺寸、网格间距、差分阶数、CFL 系数、时间步长、PML 层数、震源参数、各实验快照时刻、优化方法与系数等全部参数，用于论文复现。|
| `cfl_stability_table.png` | **CFL 稳定性表**。列出各方案的稳定时间步长上限、实际使用的时间步长、有效 CFL 数和安全裕度，确认所有模拟方案均满足稳定性条件。|

### 频散分析

| 文件名 | 含义 |
|--------|------|
| `dispersion_curves.png` | **频散误差对比图（核心图）**。上半图：数值波数比 k_num/k 随归一化波数 θ/π 的变化，理想值为 1。下半图（对数轴）：绝对频散误差 \|k_num·Δ - k·Δ\| 随波数的变化。红色为 Taylor 系数，蓝色为优化系数，黑色虚线为目标阈值 10⁻⁴。**关键结论**：优化系数的误差曲线应完全位于阈值线以下，呈现等纹波特征。|
| `sa_convergence.png` | **模拟退火收敛曲线**。上半图：温度随迭代步数的衰减过程，虚线标记 restart 重启点。下半图：目标函数值随迭代的收敛过程，红色为历史最优值，浅蓝色为当前解，绿色五角星为局部精修阶段，黑色虚线为目标阈值。用于验证优化过程收敛充分。|

### 实验 1 输出

| 文件名 | 含义 |
|--------|------|
| `exp1_regular_vs_staggered.png` | **常规网格与交错网格波场并排对比**。每行一个时刻，左列常规网格，右列交错网格。色标统一。可观察到交错网格的波前更锐利、数值频散更小。该图用于定性展示两种离散方案的差异，不作严格同阶精度比较。|

### 实验 2 输出

| 文件名 | 含义 |
|--------|------|
| `exp2_noabsorb_vs_pml.png` | **无吸收边界与 PML 边界波场快照对比**。左图为无吸收边界，右图为 PML 边界。在 t = 1.65 s 时，无吸收边界的四周可见明显的边界反射波（虚假回波），而 PML 边界的波场干净，出射波被有效吸收。|
| `exp2_energy_noabsorb_vs_pml.png` | **物理域内剩余能量对比**。左图：归一化速度平方和随时间变化。无吸收边界（蓝色）的能量在波到达边界后基本不衰减（反射回来），PML 边界（红色）的能量持续下降。右图：边缘区域能量占比随时间变化。注意：该能量指标为 vx²+vz² 速度平方和，非真正弹性总能量。|
| `exp2_pml_boundary_energy_time.png` | **PML 边界带能量随时间变化**。仅关注模型边缘窄带区域的归一化速度能量，对比无吸收与 PML 方案。PML 方案的边界带能量在波到达后快速衰减。|

### 实验 3 输出

| 文件名 | 含义 |
|--------|------|
| `exp3_standard_vs_minimax.png` | **标准系数与优化系数波场并排对比**。左列 Taylor 系数，右列优化系数。视觉上两者差异通常不大（因为整体波场幅度远大于频散差异），需要配合定量指标判断。此图仅用于定性展示。|
| `exp3_taylor_optimized_difference.png` | **波场三联对比图（含差值场）**。每行三列：左为 Taylor 系数波场，中为优化系数波场，右为差值场（优化 - Taylor）。差值场使用独立色标，显示了两组系数波场的细微差异分布。标注有相对 L2 误差值。|
| `exp3_taylor_optimized_difference_zoom.png` | **差异最大区域局部放大图**。自动定位到差值场最大的区域进行放大，从左到右为 Taylor 波场、优化波场、差值场的局部细节。用于观察优化系数在局部的改善效果。|

### 数据文件说明

| 文件名 | 含义 |
|--------|------|
| `coefficients.txt` | 优化方法、目标误差、Taylor 系数与优化系数的逐项值 |
| `summary.txt` | 全部实验结果的文本摘要报告 |
| `summary.mat` | 全部结果的 MATLAB 结构体存档 |
| `reproducibility_manifest.txt` | 可复现清单：运行模式、参数、系数、MATLAB 版本、耗时等 |
| `exp*_*.mat` | 各实验的完整波场数据（快照 + 能量曲线） |
| `*.csv` | 各类数值数据的文本导出，方便在 Excel/Python 中二次分析 |

## 技术说明

### 差分系数优化

本文所称"基于最大范数目标函数的优化系数"，是以最大绝对频散误差为目标函数：

```
E(c) = max |k_num·Δ - k·Δ|,  kh ∈ [0, 0.60π]
```

目标阈值 10⁻⁴ 与 Zhang & Yao (2013) Eq. 23 一致。优化方法为带 restart 和局部精修的模拟退火算法。

### PML 吸收边界

PML 层通过在物理模型外侧扩展网格实现（外扩方式），波场输出时裁剪回原始物理域。震源始终位于物理模型内部。

### 波场可视化色标

波场快照使用蓝-白-红发散色图（`ewm_wavefield_colormap.m`），蓝色表示负振幅，白色表示零，红色表示正振幅。色标范围使用 98.5 分位数截断，避免震源极值压制传播波前的显示对比度。

## 源代码分类索引

### 求解器（3 个）
- `ewm_staggered_solver.m` — 交错网格弹性波有限差分求解器（支持 PML）
- `ewm_regular_solver.m` — 常规网格求解器
- `ewm_pml2d.m` — PML 阻尼剖面计算

### 数值方法（5 个）
- `ewm_fd_coefficients.m` — Taylor 展开差分系数
- `ewm_optimize_minimax_coeffs.m` — 模拟退火优化
- `ewm_dispersion_metrics.m` — 频散误差评价
- `ewm_stable_dt.m` — CFL 稳定时间步长

### 模型与配置（5 个）
- `ewm_default_config.m` — 两种运行模式统一配置
- `ewm_load_marmousi.m` — Marmousi 模型加载
- `ewm_prepare_elastic_params.m` — 弹性参数计算
- `ewm_build_marmousi_10m_cache.m` — 细网格缓存构建

### 绘图（20+ 个）
- `ewm_plot_model.m` — 模型可视化
- `ewm_plot_ricker_wavelet.m` — 震源波形
- `ewm_plot_dispersion.m` — 频散曲线
- `ewm_plot_dispersion_signed.m` — 有符号频散误差
- `ewm_plot_sa_convergence.m` — SA 收敛曲线
- `ewm_plot_comparison.m` — 通用并排对比
- `ewm_plot_wavefield_triptych.m` — 三联图（含差值场）
- `ewm_plot_wavefield_zoom.m` — 局部放大图
- `ewm_plot_pml_snapshot_pair.m` — PML 快照对
- `ewm_plot_energy.m` — 能量演化曲线
- `ewm_plot_pml_boundary_energy.m` — PML 边界带能量
- `ewm_plot_exp3_metric_curves.m` — 时变误差曲线
- `ewm_plot_reflection_difference.m` — 反射残差
- `ewm_plot_pml_reflection_comparison.m` — PML 反射对比

### 工具（7 个）
- `ewm_wavefield_colormap.m` — 蓝白红发散色图
- `ewm_apply_chinese_style.m` — 中文字体设置
- `ewm_save_figure.m` — 图形保存（300 DPI）
- `ewm_pad2d.m` — 二维边界填充
- `ewm_pick_snapshot_indices.m` — 快照时间索引
- `ewm_boundary_ratio.m` — 边界能量占比
- `ewm_ricker.m` — Ricker 子波

### 输出（5 个）
- `ewm_write_summary.m` — 结果文本摘要
- `ewm_write_required_tables.m` — 参数表格（CSV + PNG）
- `ewm_light_result.m` — 轻量结果结构
- `ewm_finalize_saved_results.m` — 结果存档
- `ewm_solver_display_name.m` — 求解器命名
