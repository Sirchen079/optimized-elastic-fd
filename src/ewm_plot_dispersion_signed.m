function ewm_plot_dispersion_signed(standardCoeff, optimizedCoeff, khMax, targetError, outFile)
%EWM_PLOT_DISPERSION_SIGNED 仿师兄风格的一阶导数频散误差图。
%
% 纵轴为有符号绝对误差 k_num*Δ − kh（线性轴），横轴为奈奎斯特波数百分比
% (0–100 %)，完整展示通带内的等纹波特征和通带外的急剧衰减。
%
% 输入：
%   standardCoeff   — Taylor 交错网格系数向量（1×M）
%   optimizedCoeff  — 优化系数向量（1×M）
%   khMax           — 优化目标带宽（rad），用于画截止线
%   targetError     — 误差阈值，用于画 ±ε 参考线
%   outFile         — 输出 PNG 路径

if nargin < 3 || isempty(khMax),      khMax = 0.60 * pi;  end
if nargin < 4 || isempty(targetError), targetError = 1e-4; end
if nargin < 5 || isempty(outFile)
    outFile = 'dispersion_signed.png';
end

% --- 计算全 Nyquist 范围内的有符号误差 ---
samples = 2000;
kh = linspace(pi / samples, pi, samples).';          % kh in (0, π]
xPct = kh / pi * 100;                                % 0–100 %

stdError  = signed_error(standardCoeff,  kh);        % knum*Δ − kh
optError  = signed_error(optimizedCoeff, kh);

% --- 画图 ---
ewm_apply_chinese_style();
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 900, 560]);
ax = axes(fig);
hold(ax, 'on');

% ① 精确解：y = 0
plot(ax, [0, 100], [0, 0], 'k:', 'LineWidth', 1.2, 'DisplayName', 'Exact');

% ② Taylor 系数
plot(ax, xPct, stdError,  '-',  'Color', [0.88, 0.25, 0.05], 'LineWidth', 1.8, ...
    'DisplayName', sprintf('Conventional %d th-order staggered', 2*numel(standardCoeff)));

% ③ 优化系数
plot(ax, xPct, optError,  '-',  'Color', [0.00, 0.36, 0.90], 'LineWidth', 1.8, ...
    'DisplayName', sprintf('Optimized %d th-order staggered (SA)', 2*numel(optimizedCoeff)));

% ④ ±ε 阈值线（虚线）
yline(ax,  targetError, 'k--', 'LineWidth', 0.9, 'Alpha', 0.7, ...
    'Label', sprintf('+ε = %.0e', targetError), ...
    'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top');
yline(ax, -targetError, 'k--', 'LineWidth', 0.9, 'Alpha', 0.7, ...
    'Label', sprintf('−ε = %.0e', targetError), ...
    'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');

% ⑤ 优化带宽截止线
xcut = khMax / pi * 100;
xline(ax, xcut, ':', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.0, 'Alpha', 0.8, ...
    'Label', sprintf('k_c = %.0f%%', xcut), ...
    'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top');

% --- 坐标轴设置 ---
% 纵轴量级：先看一眼数据范围，再自动对齐到师兄那种 ×10^-3 感觉
allErr = [stdError(:); optError(:)];
yMax = max(abs(allErr(isfinite(allErr))));
yLim = min(yMax * 1.15, 3e-3);
ylim(ax, [-yLim, yLim * 0.6]);
xlim(ax, [0, 100]);

xlabel(ax, 'Percentage of Nyquist wavenumber (%)');
ylabel(ax, 'Absolute error (k_{num}\Delta - k\Delta)');
title(ax, sprintf('first-order staggered grid error  (M=%d, \\epsilon=%.0e)', ...
    numel(optimizedCoeff), targetError), 'FontSize', 13);
legend(ax, 'Location', 'southwest', 'FontSize', 10);
grid(ax, 'on');
box(ax, 'on');
set(ax, 'LineWidth', 0.9, 'FontSize', 11);

ewm_save_figure(fig, outFile);
close(fig);
end

% --- 局部函数 ---
function err = signed_error(coeff, kh)
% 有符号绝对误差：k_num*Δ − kh
offset = (0:numel(coeff)-1) + 0.5;
knum_delta = 2 * sin(kh * offset) * coeff(:);    % k_num*Δ
err = knum_delta - kh;
end
