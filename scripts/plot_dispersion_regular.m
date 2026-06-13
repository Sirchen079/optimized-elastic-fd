%PLOT_DISPERSION_REGULAR 独立脚本：仿师兄那张图的版式，画 10 阶交错网格
% 一阶差分误差曲线（中文标注）。
%
% 横轴：奈奎斯特波数百分比 (0–100 %)
% 纵轴：有符号绝对误差  k_num*Δ − kh
%
% 系数来源：results_standard/coefficients.txt
%   - Taylor   10 阶交错系数（常规）
%   - 优化     10 阶交错系数（模拟退火 minimax）
%
% 用法：MATLAB 中直接运行本文件，输出图保存到
%   results_standard/figures/dispersion_staggered_zh.png

clear; clc;
projectDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectDir, 'src'));

% -------- 1. 读系数 --------
coeffFile = fullfile(projectDir, 'results_standard', 'coefficients.txt');
[taylorCoeff, optimizedCoeff] = read_coeff_file(coeffFile);

order = numel(optimizedCoeff);   % 单侧系数个数，对应 2*order 阶
fprintf('Taylor   %d 阶系数：%s\n', 2*order, sprintf('%.6g  ', taylorCoeff));
fprintf('优化     %d 阶系数：%s\n', 2*order, sprintf('%.6g  ', optimizedCoeff));

% -------- 2. 计算频散误差 --------
samples = 4000;
kh   = linspace(pi/samples, pi, samples).';
xPct = kh / pi * 100;

errTaylor    = staggered_signed_error(taylorCoeff,    kh);
errOptimized = staggered_signed_error(optimizedCoeff, kh);

% -------- 3. 画图 --------
ewm_apply_chinese_style();
fig = figure('Color', 'w', 'Position', [100, 100, 980, 680]);
ax  = axes(fig); hold(ax, 'on');

% 精确解：y = 0（黑点线）
plot(ax, [0, 100], [0, 0], 'k:', 'LineWidth', 1.4, 'DisplayName', '精确解');

% 优化曲线（红色实线）
plot(ax, xPct, errOptimized, '-', 'Color', [0.86, 0.10, 0.10], ...
     'LineWidth', 1.8, 'DisplayName', sprintf('优化 %d 阶', 2*order));

% Taylor 曲线（黑色实线）
plot(ax, xPct, errTaylor, 'k-', 'LineWidth', 1.5, ...
     'DisplayName', sprintf('常规 %d 阶', 2*order));

% ±1e-4 误差阈值线（黑色虚线，带标注）
targetErr = 1e-4;
yline(ax,  targetErr, 'k--', 'LineWidth', 1.0, ...
    'Label', sprintf('+%.0e', targetErr), ...
    'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top', ...
    'FontSize', 10, 'HandleVisibility', 'off');
yline(ax, -targetErr, 'k--', 'LineWidth', 1.0, ...
    'Label', sprintf('\x2212%.0e', targetErr), ...
    'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', ...
    'FontSize', 10, 'HandleVisibility', 'off');

xlim(ax, [0, 100]);
ylim(ax, [-2e-3, 1e-3]);
xlabel(ax, '奈奎斯特波数百分比 (%)');
ylabel(ax, '绝对误差');
title(ax, sprintf('一阶交错网格差分误差（%d 阶）', 2*order), 'FontSize', 14);
legend(ax, 'Location', 'northeast', 'FontSize', 11, 'Box', 'on');
grid(ax, 'on'); box(ax, 'on');
set(ax, 'LineWidth', 0.9, 'FontSize', 12);

% -------- 4. 保存 --------
outFile = fullfile(projectDir, 'results_standard', 'figures', ...
                   'dispersion_staggered_zh.png');
ewm_save_figure(fig, outFile);
fprintf('\n图已保存到：%s\n', outFile);
close(fig);


% ====================== 局部函数 ======================

function err = staggered_signed_error(coeff, kh)
% 交错网格一阶导数的有符号绝对误差：
%   k_num*Δ = 2 * Σ_{n=0..N-1} c_{n+1} * sin((n+0.5)*kh)
%   err     = k_num*Δ − kh
offset = (0:numel(coeff)-1) + 0.5;
knum   = 2 * sin(kh * offset) * coeff(:);
err    = knum - kh;
end

function [taylorCoeff, optimizedCoeff] = read_coeff_file(coeffFile)
taylorCoeff    = [];
optimizedCoeff = [];

if ~isfile(coeffFile)
    error('未找到系数文件：%s', coeffFile);
end

fid = fopen(coeffFile, 'r');
cleaner = onCleanup(@() fclose(fid));
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if startsWith(line, 'Taylor系数 =')
        taylorCoeff = sscanf(line(length('Taylor系数 =') + 1:end), '%f').';
    elseif startsWith(line, '优化系数 =')
        optimizedCoeff = sscanf(line(length('优化系数 =') + 1:end), '%f').';
    end
end

if isempty(taylorCoeff) || isempty(optimizedCoeff)
    error('coefficients.txt 中未找到 "Taylor系数 =" 或 "优化系数 =" 行。');
end
end
