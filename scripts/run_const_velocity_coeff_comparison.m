function results = run_const_velocity_coeff_comparison()
%RUN_CONST_VELOCITY_COEFF_COMPARISON 常速度模型下 Taylor 系数与
%基于最大范数目标函数的优化系数的单道波形对比。
%
% - 模型：均匀弹性介质（Vp、Vs、密度处处相同），6 km × 6 km 网格。
% - 边界：四周 PML 吸收，避免人工反射干扰单道波形。
% - 源：Ricker 子波垂直点力（Fz）置于模型中央。
%   设计意图：力方向 z̄ 与传播方向 x̄ 垂直，P 辐射系数 cos(90°)=0、
%   S 辐射系数 sin(90°)=1，水平偏移接收点上 vz 道得到纯 S 直达，
%   主频对应 kh_S=0.5π 恰好落在 Taylor 误差陡升区（~3×10⁻³），
%   而优化系数仍 <10⁻⁴，累积 2 km（20 个 S 波长）相位差异肉眼可辨。
% - 接收：源右侧 2 km 处提取 vz 单道时间序列。
% - 输出：全时段单道叠加图（Taylor 蓝实线 / 优化 红虚线）。

projectDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectDir, 'src'));
ewm_apply_chinese_style();

outDir = fullfile(projectDir, 'results_const_velocity');
figDir = fullfile(outDir, 'figures');
ewm_ensure_dir(outDir);
ewm_ensure_dir(figDir);

% ------------- 1. 构造均匀常速度弹性模型 -------------
% 选用 Vp/Vs = 1.5、dx = 25 m、f0 = 20 Hz：
%   - 峰值 kh_p = 2π·20·25/3000 ≈ 0.33π，P 波 Taylor 与优化都能良好支持
%   - 峰值 kh_s = 2π·20·25/2000 ≈ 0.50π，落在 Taylor 5 阶系数误差迅速
%     上升（~3×10⁻³）的区段，而优化系数仍 <10⁻⁴
%   - 累积传播 1.50 s + 2.00 km 偏移让单道差异充分显化
% 域扩到 8 km × 8 km：源/收点居中后，接收点距右侧 PML 仍有 2 km，
% 把 PML 反射延迟推到 (2·2km + 2km)/Vs ≈ 3 s，远远超出 1.5 s 模拟窗口，
% 单道里看到的差异就只来自频散与源耦合。
model = build_constant_model( ...
    'vp', 3000, 'vs', 2000, 'rho', 2300, ...
    'dx', 25, 'nx', 321, 'nz', 321);
fprintf('常速度模型：Vp = %.0f m/s, Vs = %.0f m/s, ρ = %.0f kg/m³\n', ...
    3000, 2000, 2300);
fprintf('  网格：%d × %d，dx = dz = %.1f m，区域 %.1f km × %.1f km\n', ...
    model.nz, model.nx, model.dx, (model.nx-1)*model.dx/1000, ...
    (model.nz-1)*model.dz/1000);

% ------------- 2. 差分系数 -------------
order = 5;
standardCoeff = ewm_fd_coefficients('standard', order);
% 论文里基于最大范数目标函数（Zhang & Yao 2013 Eq.23）的优化系数，
% 来自 results_standard/coefficients.txt：直接复用以避免重新跑 SA。
optimizedCoeff = [ ...
    1.23717720156044, ...
   -0.108727049806825, ...
    0.0237950229290865, ...
   -0.0052569485184926, ...
    0.000758608093783266];

fprintf('Taylor 系数：%s\n', sprintf('%+.6f ', standardCoeff));
fprintf('优化系数  ：%s\n', sprintf('%+.6f ', optimizedCoeff));

% ------------- 3. 仿真参数 -------------
cfl = 0.42;
dtT = ewm_stable_dt(model, standardCoeff, cfl);
dtO = ewm_stable_dt(model, optimizedCoeff, cfl);
dt = min(dtT, dtO);
totalTime = 1.50;          % 秒
nt = ceil(totalTime / dt) + 1;

% 源、接收
sourceZkm = (model.nz - 1) * model.dz / 1000 / 2;   % 模型中心深度（3.0 km）
sourceXkm = (model.nx - 1) * model.dx / 1000 / 2;   % 模型中心水平（3.0 km）
recOffsetKm = 2.00;                                  % 偏移 2 km，足够 S 波频散累积
recZkm = sourceZkm;
recXkm = sourceXkm + recOffsetKm;

recIz = round(recZkm * 1000 / model.dz) + 1;
recIx = round(recXkm * 1000 / model.dx) + 1;

sim = struct();
sim.dt = dt;
sim.nt = nt;
sim.nPml = 20;
sim.f0 = 20.0;
sim.sourceDelayCycles = 1.5;
sim.sourceAmplitude = 1.0;
sim.sourceType = 'pointForceZ';   % 垂直点力，水平偏移接收点上得到纯 S 直达
sim.sourceDepthM = sourceZkm * 1000;
sim.sourceXFraction = sourceXkm / ((model.nx - 1) * model.dx / 1000);

% 空间高斯 stamp：σ = 1 个网格点的 5×5 加权核，sum=1。
% 作用：把 δ 点源在网格上注入的 kh > 0.7π 短波长成分压下去，
% 同时保留 kh ≤ 0.5π 的有效信号（主频 S 波 λ = 100 m = 4·dx）。
sim.sourceStamp = build_gaussian_stamp(2, 1.0);

fprintf('震源类型：垂直点力 Fz（高斯空间平滑，σ=1 网格点，footprint 5×5）\n');
sim.snapshotFractions = [0.35, 0.65, 0.95];
sim.energyStride = max(1, round(nt / 200));
sim.outputDir = [];
sim.receivers.iz = recIz;
sim.receivers.ix = recIx;

fprintf('源位置：(z=%.2f km, x=%.2f km)，接收点：(z=%.2f km, x=%.2f km)，偏移 %.2f km\n', ...
    sourceZkm, sourceXkm, recZkm, recXkm, recOffsetKm);
fprintf('dt = %.3e s，nt = %d（总时长 %.3f s），f0 = %.1f Hz\n', ...
    dt, nt, (nt-1)*dt, sim.f0);

% ------------- 4. 跑两次正演 -------------
fprintf('\n[1/2] 跑 Taylor 系数...\n');
resultTaylor = ewm_staggered_solver(model, sim, standardCoeff, true, ...
    'const_velocity_taylor');
fprintf('[2/2] 跑优化系数...\n');
resultOpt = ewm_staggered_solver(model, sim, optimizedCoeff, true, ...
    'const_velocity_optimized');

save(fullfile(outDir, 'const_velocity_taylor.mat'), 'resultTaylor', '-v7.3');
save(fullfile(outDir, 'const_velocity_optimized.mat'), 'resultOpt', '-v7.3');

% ------------- 5. 出图：单道时间序列叠加 -------------
traceFig = fullfile(figDir, 'const_velocity_trace_compare.png');
plot_trace_compare(resultTaylor, resultOpt, ...
    'Taylor 系数', '基于最大范数目标函数的优化系数', ...
    [recZkm, recXkm], traceFig);

% ------------- 6. 保存单道数据 -------------
traceCsv = fullfile(outDir, 'const_velocity_trace.csv');
write_trace_csv(traceCsv, resultTaylor, resultOpt);

% ------------- 7. 导出实验参数 markdown -------------
mdFile = fullfile(outDir, 'const_velocity_experiment_parameters.md');
write_parameters_markdown(mdFile, model, sim, cfl, dtT, dtO, ...
    standardCoeff, optimizedCoeff, order, ...
    [sourceZkm, sourceXkm], [recZkm, recXkm, recOffsetKm], ...
    {traceFig, traceCsv});

results = struct();
results.model = struct('vp', 3000, 'vs', 2000, 'rho', 2300, ...
    'dx', model.dx, 'dz', model.dz, 'nx', model.nx, 'nz', model.nz);
results.sim = sim;
results.standardCoeff = standardCoeff;
results.optimizedCoeff = optimizedCoeff;
results.receiver = struct('zKm', recZkm, 'xKm', recXkm, ...
    'iz', recIz, 'ix', recIx, 'offsetKm', recOffsetKm);
results.figureFiles = struct('trace', traceFig);
save(fullfile(outDir, 'summary.mat'), 'results', '-v7.3');

fprintf('\n输出文件：\n');
fprintf('  单道对比图：%s\n', traceFig);
fprintf('  单道 CSV  ：%s\n', traceCsv);
fprintf('  参数 MD   ：%s\n', mdFile);
fprintf('  汇总 mat  ：%s\n', fullfile(outDir, 'summary.mat'));
end

% ============================================================
function stamp = build_gaussian_stamp(halfWidth, sigmaGridPts)
%构造空间高斯权重核，footprint 为 (2·halfWidth+1) × (2·halfWidth+1)，
%权重和归一化为 1。sigmaGridPts 单位为网格点。
[I, J] = ndgrid(-halfWidth:halfWidth, -halfWidth:halfWidth);
stamp = exp(-(I.^2 + J.^2) / (2 * sigmaGridPts^2));
stamp = stamp / sum(stamp(:));
end

% ============================================================
function model = build_constant_model(varargin)
%构造均匀常速度弹性模型
p = inputParser;
addParameter(p, 'vp', 3000);
addParameter(p, 'vs', 1800);
addParameter(p, 'rho', 2300);
addParameter(p, 'dx', 20);
addParameter(p, 'nx', 251);
addParameter(p, 'nz', 251);
parse(p, varargin{:});
r = p.Results;

model = struct();
model.vp = r.vp * ones(r.nz, r.nx);
model.vs = r.vs * ones(r.nz, r.nx);
model.rho = r.rho * ones(r.nz, r.nx);
model.dx = r.dx;
model.dz = r.dx;
model.nx = r.nx;
model.nz = r.nz;
model.x = (0:r.nx-1) * r.dx;
model.z = (0:r.nz-1) * r.dx;
model.source = '常速度均匀模型（自构造）';
end

% ============================================================
function plot_trace_compare(leftResult, rightResult, ...
    leftLabel, rightLabel, recZX, outFile)
%绘制单道全时段时间序列叠加对比
ewm_apply_chinese_style();
t = leftResult.trace.time;
vzL = leftResult.trace.vz(:, 1);
vzR = rightResult.trace.vz(:, 1);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1500, 560]);
ax = axes(fig);
set(ax, 'FontSize', 22);
hold(ax, 'on'); grid(ax, 'on');
plot(ax, t, vzL, 'b-',  'LineWidth', 1.8, 'DisplayName', leftLabel);
plot(ax, t, vzR, 'r--', 'LineWidth', 1.6, 'DisplayName', rightLabel);
xlabel(ax, '时间 (s)', 'FontSize', 24);
ylabel(ax, '质点垂直振动速度 v_z (m·s^{-1})', 'FontSize', 24);
title(ax, sprintf('接收点 (z = %.2f km, x = %.2f km) 单道时间序列对比', ...
    recZX(1), recZX(2)), 'FontSize', 24);
legend(ax, 'Location', 'best', 'FontSize', 20);
xlim(ax, [t(1), t(end)]);

ewm_save_figure(fig, outFile);
close(fig);
end

% ============================================================
function write_parameters_markdown(outFile, model, sim, cfl, dtT, dtO, ...
    standardCoeff, optimizedCoeff, order, srcZX, recZXOffset, outputPaths)
%将本次常速度对比实验的全部参数写入 markdown 文件。
fid = fopen(outFile, 'w', 'n', 'UTF-8');
cleanup = onCleanup(@() fclose(fid));

vp = model.vp(1, 1);
vs = model.vs(1, 1);
rho = model.rho(1, 1);
dx = model.dx; dz = model.dz;
nx = model.nx; nz = model.nz;
xKm = (nx - 1) * dx / 1000;
zKm = (nz - 1) * dz / 1000;
nt = sim.nt; dt = sim.dt;
f0 = sim.f0; nPml = sim.nPml;
sourceZkm = srcZX(1); sourceXkm = srcZX(2);
recZkm = recZXOffset(1); recXkm = recZXOffset(2); recOffsetKm = recZXOffset(3);
totalTime = (nt - 1) * dt;

% 计算频散关键指标
khP = 2 * pi * f0 * dx / vp;
khS = 2 * pi * f0 * dx / vs;
stencilGainT = sum(abs(standardCoeff));
stencilGainO = sum(abs(optimizedCoeff));

fprintf(fid, '# 常速度模型差分系数对比实验参数\n\n');
fprintf(fid, '本文件由 `run_const_velocity_coeff_comparison.m` 自动生成。\n\n');

fprintf(fid, '## 1. 介质模型（均匀常速度弹性体）\n\n');
fprintf(fid, '| 参数 | 数值 | 说明 |\n');
fprintf(fid, '| --- | --- | --- |\n');
fprintf(fid, '| 纵波速度 Vp | %.0f m/s | 全网格均匀 |\n', vp);
fprintf(fid, '| 横波速度 Vs | %.0f m/s | Vp/Vs = %.3f |\n', vs, vp / vs);
fprintf(fid, '| 密度 ρ | %.0f kg/m³ | 全网格均匀 |\n', rho);
fprintf(fid, '| 网格间距 dx = dz | %.1f m | 各向同性 |\n', dx);
fprintf(fid, '| 网格点数 nx × nz | %d × %d | — |\n', nx, nz);
fprintf(fid, '| 物理区域 | %.2f km × %.2f km | 水平 × 深度 |\n\n', xKm, zKm);

fprintf(fid, '## 2. 差分系数（阶数 M = %d，交错网格一阶导）\n\n', order);
fprintf(fid, '| 系数 | Taylor | 基于最大范数目标函数的优化 |\n');
fprintf(fid, '| --- | --- | --- |\n');
for k = 1:numel(standardCoeff)
    fprintf(fid, '| c_%d | %+.12g | %+.12g |\n', ...
        k, standardCoeff(k), optimizedCoeff(k));
end
fprintf(fid, '| Σ\\|c_n\\|（模板增益） | %.6f | %.6f |\n\n', ...
    stencilGainT, stencilGainO);
fprintf(fid, '优化系数取自 `results_standard/coefficients.txt`，khMax = 0.6π，目标误差 1×10⁻⁴。\n\n');

fprintf(fid, '## 3. 时间步与采样\n\n');
fprintf(fid, '| 参数 | 数值 | 说明 |\n');
fprintf(fid, '| --- | --- | --- |\n');
fprintf(fid, '| CFL 数 | %.3f | 经验稳定上限 |\n', cfl);
fprintf(fid, '| Taylor 系数稳定 dt 上限 | %.6e s | — |\n', dtT);
fprintf(fid, '| 优化系数稳定 dt 上限 | %.6e s | — |\n', dtO);
fprintf(fid, '| 实际采用 dt | %.6e s | 取两者最小值 |\n', dt);
fprintf(fid, '| 总时间步数 nt | %d | — |\n', nt);
fprintf(fid, '| 模拟总时长 | %.4f s | (nt−1)·dt |\n\n', totalTime);

fprintf(fid, '## 4. 震源\n\n');
fprintf(fid, '| 参数 | 数值 | 说明 |\n');
fprintf(fid, '| --- | --- | --- |\n');
fprintf(fid, '| 子波类型 | Ricker | — |\n');
fprintf(fid, '| 加载方式 | 垂直点力 Fz | 力方向 z̄ 与传播方向 x̄ 垂直：P 辐射 ∝ cos(90°)=0，S 辐射 ∝ sin(90°)=1 |\n');
fprintf(fid, '| 空间分布 | 高斯 stamp σ=1·dx，5×5 网格 | 抑制 δ 源在 kh>0.7π 区注入的不稳定短波长 |\n');
fprintf(fid, '| 主频 f₀ | %.2f Hz | — |\n', f0);
fprintf(fid, '| 延迟周期数 | %.2f | t₀ = %.4f s |\n', ...
    sim.sourceDelayCycles, sim.sourceDelayCycles / f0);
fprintf(fid, '| 振幅 | %.3g | — |\n', sim.sourceAmplitude);
fprintf(fid, '| 位置 (z, x) | (%.3f km, %.3f km) | 模型中心 |\n\n', ...
    sourceZkm, sourceXkm);

fprintf(fid, '## 5. 接收点\n\n');
fprintf(fid, '| 参数 | 数值 | 说明 |\n');
fprintf(fid, '| --- | --- | --- |\n');
fprintf(fid, '| 位置 (z, x) | (%.3f km, %.3f km) | 源右侧 |\n', recZkm, recXkm);
fprintf(fid, '| 距源偏移 | %.3f km | 水平直达 |\n', recOffsetKm);
fprintf(fid, '| 累积 S 波长数 | %.1f | offset · f₀ / Vs |\n', recOffsetKm * 1000 * f0 / vs);
fprintf(fid, '| 记录分量 | vz（主分析）、vx（备查） | 全时间步采样 |\n');
fprintf(fid, '| 单道物理意义 | 纯 S 直达 | 垂直力 + 水平传播几何下 vz 等于 S 横向位移分量 |\n\n');

fprintf(fid, '## 6. 边界条件\n\n');
pmlBufferKm = (xKm - recXkm) - nPml * dx / 1000;
pmlEchoTime = (2 * (xKm - recXkm) * 1000) / vs + sim.sourceDelayCycles / f0;
fprintf(fid, '- 四周 PML 吸收，层数 = %d。\n', nPml);
fprintf(fid, '- 接收点距右侧 PML 边缘缓冲 = %.2f km。\n', pmlBufferKm);
fprintf(fid, '- 右侧 PML 反射回到接收点的最早到时 ≈ %.2f s，远大于总模拟时长 %.2f s；窗口内不会被污染。\n\n', ...
    pmlEchoTime, totalTime);

fprintf(fid, '## 7. 关键频散指标\n\n');
fprintf(fid, '主频 f₀ 处的归一化波数 kh（k·Δx）：\n\n');
fprintf(fid, '| 波类 | kh = 2πf₀·Δx/v | 占 π 的比例 |\n');
fprintf(fid, '| --- | --- | --- |\n');
fprintf(fid, '| P 波 | %.4f | %.3fπ |\n', khP, khP / pi);
fprintf(fid, '| S 波 | %.4f | %.3fπ |\n\n', khS, khS / pi);
fprintf(fid, '> 垂直点力 + 水平偏移接收：vz 道为纯 S 直达。S 波主频 kh = 0.5π 落在 Taylor\n');
fprintf(fid, '> 5 阶系数相位误差陡升区（~3×10⁻³），而优化系数 <10⁻⁴；2 km 偏移累积约 20\n');
fprintf(fid, '> 个 S 波长，Taylor 累积相位误差约 0.4 rad，单道尾部应显出明显的高频频散环。\n\n');

fprintf(fid, '## 8. 输出文件\n\n');
for k = 1:numel(outputPaths)
    [~, name, ext] = fileparts(outputPaths{k});
    fprintf(fid, '- `%s%s`\n', name, ext);
end
end

% ============================================================
function write_trace_csv(outFile, leftResult, rightResult)
fid = fopen(outFile, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'time_s,vz_taylor,vz_optimized,vx_taylor,vx_optimized\n');
t = leftResult.trace.time;
vzL = leftResult.trace.vz(:, 1);
vzR = rightResult.trace.vz(:, 1);
vxL = leftResult.trace.vx(:, 1);
vxR = rightResult.trace.vx(:, 1);
for k = 1:numel(t)
    fprintf(fid, '%.9g,%.9g,%.9g,%.9g,%.9g\n', ...
        t(k), vzL(k), vzR(k), vxL(k), vxR(k));
end
end
