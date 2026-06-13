%PLOT_DISPERSION_HIGH_WAVENUMBER 独立脚本：高波数域优化系数与泰勒系数误差对比。
%
% 横轴：归一化波数 kh/pi (高波数域 0.5–1.0)
% 纵轴：有符号绝对误差 k_num*Delta - kh
%
% 系数来源：results_standard/coefficients.txt
%   - Taylor   10 阶交错系数 (常规)
%   - 优化     10 阶交错系数 (模拟退火 minimax)
%
% 输出：results_standard/figures/dispersion_high_wavenumber.png
%
% 用法：MATLAB 中直接运行本文件。

clear; clc;
projectDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectDir, 'src'));

% -------- 1. 读系数 --------
coeffFile = fullfile(projectDir, 'results_standard', 'coefficients.txt');
[taylorCoeff, optimizedCoeff, khMax, targetError] = read_coeff_file(coeffFile);

order = numel(optimizedCoeff);
fprintf('Taylor   %d 阶系数：%s\n', 2*order, sprintf('%.6g  ', taylorCoeff));
fprintf('优化     %d 阶系数：%s\n', 2*order, sprintf('%.6g  ', optimizedCoeff));
fprintf('优化带宽 khMax = %.4f (%.1f%% Nyquist)\n', khMax, khMax/pi*100);

% -------- 2. 计算高波数域误差 --------
samples = 4000;
kh = linspace(pi/samples, pi, samples).';
xNorm = kh / pi;  % 0–1

errTaylor    = staggered_signed_error(taylorCoeff,    kh);
errOptimized = staggered_signed_error(optimizedCoeff, kh);

% 高波数域统计 (kh/pi >= 0.5)
hwMask = xNorm >= 0.5;
khHw = kh(hwMask);
errTaylorHw = errTaylor(hwMask);
errOptimizedHw = errOptimized(hwMask);

maxAbsTaylorHw    = max(abs(errTaylorHw));
maxAbsOptimizedHw = max(abs(errOptimizedHw));
rmsTaylorHw       = rms(errTaylorHw);
rmsOptimizedHw    = rms(errOptimizedHw);

fprintf('\n--- 高波数域 (kh/pi >= 0.5) 误差统计 ---\n');
fprintf('Taylor   最大绝对误差 = %.6g\n', maxAbsTaylorHw);
fprintf('优化     最大绝对误差 = %.6g\n', maxAbsOptimizedHw);
fprintf('Taylor   RMS 误差     = %.6g\n', rmsTaylorHw);
fprintf('优化     RMS 误差     = %.6g\n', rmsOptimizedHw);
fprintf('最大误差改善倍数       = %.2f\n', maxAbsTaylorHw / maxAbsOptimizedHw);
fprintf('RMS 误差改善倍数       = %.2f\n', rmsTaylorHw / rmsOptimizedHw);

% -------- 3. 画图 --------
ewm_apply_chinese_style();
fig = figure('Color', 'w', 'Position', [100, 100, 1000, 820]);

xCut = khMax / pi;

% 用 tiledlayout 控制间距，上方留空给总标题
tlo = tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
tlo.Title.String = sprintf(['高波数域优化系数与泰勒系数误差对比 | ', ...
    '高波数域 (kh/\\pi \\geq 0.5): Taylor max|err| = %.4g, ', ...
    '优化 max|err| = %.4g, 改善 %.1f 倍'], ...
    maxAbsTaylorHw, maxAbsOptimizedHw, maxAbsTaylorHw / maxAbsOptimizedHw);
tlo.Title.FontWeight = 'bold';
tlo.Title.FontSize = 13;
tlo.Title.Interpreter = 'tex';

% ---- 上子图：有符号绝对误差 (高波数域聚焦) ----
ax1 = nexttile(1); hold(ax1, 'on');

plot(ax1, [0.5, 1.0], [0, 0], 'k:', 'LineWidth', 1.2, 'DisplayName', '精确解');

plot(ax1, xNorm, errTaylor, '-', 'Color', [0.88, 0.25, 0.05], 'LineWidth', 1.6, ...
    'DisplayName', sprintf('常规 %d 阶 (Taylor)', 2*order));

plot(ax1, xNorm, errOptimized, '-', 'Color', [0.00, 0.36, 0.90], 'LineWidth', 2.0, ...
    'DisplayName', sprintf('优化 %d 阶 (SA minimax)', 2*order));

xline(ax1, xCut, '--', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.2, ...
    'DisplayName', sprintf('优化带宽 k_c = %.2f\\pi', xCut));

xline(ax1, 0.5, ':', 'Color', [0.6, 0.2, 0.6], 'LineWidth', 1.0, ...
    'DisplayName', '高波数域起点 0.5\\pi');

xlim(ax1, [0.5, 1.0]);
hwErrAll = [errTaylorHw(:); errOptimizedHw(:)];
yMaxHw = max(abs(hwErrAll(isfinite(hwErrAll))));
ylim(ax1, [-yMaxHw * 1.1, yMaxHw * 1.1]);
xlabel(ax1, '归一化波数 kh / \\pi');
ylabel(ax1, '绝对误差 k_{num}\\Delta - k\\Delta');
title(ax1, sprintf('高波数域频散误差对比 (%d 阶交错网格)', 2*order), 'FontSize', 12);
legend(ax1, 'Location', 'southwest', 'FontSize', 10, 'Box', 'on');
grid(ax1, 'on'); box(ax1, 'on');
set(ax1, 'LineWidth', 0.9, 'FontSize', 11);

% ---- 下子图：相对误差 |k_num/k - 1| (对数轴，高波数域) ----
ax2 = nexttile(2); hold(ax2, 'on');

relErrTaylor    = abs(errTaylor ./ kh);
relErrOptimized = abs(errOptimized ./ kh);

floorRel = 1e-6;
relErrTaylorPlot    = max(relErrTaylor,    floorRel);
relErrOptimizedPlot = max(relErrOptimized, floorRel);

semilogy(ax2, xNorm, relErrTaylorPlot, '-', 'Color', [0.88, 0.25, 0.05], 'LineWidth', 1.6, ...
    'DisplayName', sprintf('常规 %d 阶 (Taylor)', 2*order));
semilogy(ax2, xNorm, relErrOptimizedPlot, '-', 'Color', [0.00, 0.36, 0.90], 'LineWidth', 2.0, ...
    'DisplayName', sprintf('优化 %d 阶 (SA minimax)', 2*order));

xline(ax2, xCut, '--', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.2);
xline(ax2, 0.5, ':', 'Color', [0.6, 0.2, 0.6], 'LineWidth', 1.0);

xlim(ax2, [0.5, 1.0]);
xlabel(ax2, '归一化波数 kh / \\pi');
ylabel(ax2, '相对误差 |k_{num}/k - 1|');
title(ax2, '高波数域相对误差 (对数轴)', 'FontSize', 12);
legend(ax2, 'Location', 'southwest', 'FontSize', 10, 'Box', 'on');
grid(ax2, 'on'); box(ax2, 'on');
set(ax2, 'LineWidth', 0.9, 'FontSize', 11);

% -------- 4. 保存 --------
outDir = fullfile(projectDir, 'results_standard', 'figures');
ewm_ensure_dir(outDir);
outFile = fullfile(outDir, 'dispersion_high_wavenumber.png');
ewm_save_figure(fig, outFile);
fprintf('\nFigure saved to: %s\n', outFile);
close(fig);


% ====================== 局部函数 ======================

function err = staggered_signed_error(coeff, kh)
offset = (0:numel(coeff)-1) + 0.5;
knum   = 2 * sin(kh * offset) * coeff(:);
err    = knum - kh;
end

function [taylorCoeff, optimizedCoeff, khMax, targetError] = read_coeff_file(coeffFile)
taylorCoeff    = [];
optimizedCoeff = [];
khMax          = 0.60 * pi;
targetError    = 1e-4;

if ~isfile(coeffFile)
    error('Coefficient file not found: %s', coeffFile);
end

fid = fopen(coeffFile, 'r');
cleaner = onCleanup(@() fclose(fid));
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if startsWith(line, 'Taylor系数 =')
        taylorCoeff = sscanf(line(length('Taylor系数 =') + 1:end), '%f').';
    elseif startsWith(line, '优化系数 =')
        optimizedCoeff = sscanf(line(length('优化系数 =') + 1:end), '%f').';
    elseif startsWith(line, '最大kh =')
        khMax = sscanf(line(length('最大kh =') + 1:end), '%f');
    elseif startsWith(line, '目标误差 =')
        targetError = sscanf(line(length('目标误差 =') + 1:end), '%f');
    end
end

if isempty(taylorCoeff) || isempty(optimizedCoeff)
    error('Could not find Taylor or optimized coefficients in file.');
end
end
