%PLOT_DISPERSION_ERROR_COMPARISON 独立脚本：通带等纹波 + 全波数域误差对比。
%
% 生成两张图：
%   图A: 通带聚焦等纹波图 — 线性轴，展示 SA minimax 在通带内压制误差的等纹波特征
%   图B: 全波数域绝对误差图 — 对数轴，展示完整 Nyquist 范围的误差全景
%
% 误差标准统一为 SA 目标函数：
%   有符号绝对误差: k_num*Delta - k*Delta，目标 epsilon = 1e-4
%
% 系数来源：results_standard/coefficients.txt
% 输出：
%   results_standard/figures/dispersion_passband_equiripple.png
%   results_standard/figures/dispersion_fullband_logerror.png
%
% 用法：MATLAB 中直接运行本文件。

clear; clc;
projectDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectDir, 'src'));

% ======================================================================
% 1. 读取已有系数数据
% ======================================================================
coeffFile = fullfile(projectDir, 'results_standard', 'coefficients.txt');
[taylorCoeff, optimizedCoeff, khMax, targetError] = read_coeff_file(coeffFile);

order = numel(optimizedCoeff);  % 单侧系数个数，对应 2*order 阶
fprintf('=== 系数读取 ===\n');
fprintf('Taylor   %d 阶: %s\n', 2*order, sprintf('%.6g  ', taylorCoeff));
fprintf('优化     %d 阶: %s\n', 2*order, sprintf('%.6g  ', optimizedCoeff));
fprintf('优化带宽 khMax = %.4f (%.1f%%% Nyquist)\n', khMax, khMax/pi*100);
fprintf('目标误差 epsilon = %.0e\n', targetError);

% ======================================================================
% 2. 统一误差标准计算
%    E(c) = max | k_num(c)*Delta - k*Delta |, kh in (0, khMax]
% ======================================================================
samples = 4000;
kh = linspace(pi/samples, pi, samples).';
xNorm = kh / pi;  % 0 ~ 1

% 有符号绝对误差 (SA 目标函数形式)
errTaylor    = staggered_signed_error(taylorCoeff,    kh);
errOptimized = staggered_signed_error(optimizedCoeff, kh);

% 通带内最大绝对误差 (kh <= khMax)
inBand = kh <= khMax;
E_taylor = max(abs(errTaylor(inBand)));
E_opt    = max(abs(errOptimized(inBand)));
improvement = E_taylor / E_opt;

fprintf('\n=== 126 倍数值计算依据 ===\n');
fprintf('误差标准: E(c) = max |k_num(c)*Delta - k*Delta|, kh in (0, %.4f]\n', khMax);
fprintf('E_taylor = %.15g  (Taylor %d阶通带最大绝对误差)\n', E_taylor, 2*order);
fprintf('E_opt    = %.15g  (优化系数通带最大绝对误差)\n', E_opt);
fprintf('epsilon  = %.15g  (SA 优化目标阈值)\n', targetError);
fprintf('改善倍数 = E_taylor / E_opt = %.15g / %.15g = %.3f\n', ...
    E_taylor, E_opt, improvement);
fprintf('验证: E_opt (%.3e) < epsilon (%.0e) => %s\n', ...
    E_opt, targetError, ternary(E_opt < targetError, 'PASS', 'FAIL'));

% ======================================================================
% 3. 图A：通带聚焦等纹波图
% ======================================================================
fprintf('\n=== 生成图A：通带聚焦等纹波图 ===\n');

ewm_apply_chinese_style();
xCut = khMax / pi;   % 优化带宽截止 (归一化)
xMaxA = 0.65;        % 略超通带以展示过渡

figA = figure('Color', 'w', 'Position', [100, 100, 1000, 700]);

% 用 tiledlayout，上方留空给总标题
tloA = tiledlayout(figA, 1, 1, 'Padding', 'compact');
tloA.Title.String = sprintf(['通带聚焦：优化系数 vs 泰勒系数频散误差对比 ' ...
    '(%d阶交错网格, \\epsilon = 10^{-4})'], 2*order);
tloA.Title.FontWeight = 'bold';
tloA.Title.FontSize = 13;
tloA.Title.Interpreter = 'tex';

axA = nexttile(1); hold(axA, 'on');

% 精确解 y = 0
plot(axA, [0, xMaxA], [0, 0], 'k:', 'LineWidth', 1.2, 'DisplayName', '精确解');

% Taylor 系数曲线
plot(axA, xNorm, errTaylor, '-', 'Color', [0.88, 0.25, 0.05], 'LineWidth', 1.6, ...
    'DisplayName', sprintf('常规 %d 阶 (Taylor)', 2*order));

% 优化系数曲线
plot(axA, xNorm, errOptimized, '-', 'Color', [0.00, 0.36, 0.90], 'LineWidth', 2.2, ...
    'DisplayName', sprintf('优化 %d 阶 (SA minimax)', 2*order));

% +/- epsilon 参考线
yline(axA,  targetError, 'k--', 'LineWidth', 1.0, ...
    'DisplayName', sprintf('+\\epsilon = 10^{-4}'));
yline(axA, -targetError, 'k--', 'LineWidth', 1.0, ...
    'DisplayName', sprintf('-\\epsilon = -10^{-4}'), 'HandleVisibility', 'off');

% 优化带宽截止线
xline(axA, xCut, '--', 'Color', [0.3, 0.3, 0.3], 'LineWidth', 1.2, ...
    'DisplayName', sprintf('优化带宽 k_c = %.2f\\pi', xCut));

% --- 坐标轴 ---
xlim(axA, [0, xMaxA]);
% Y 轴范围：把 Taylor 的最大误差纳入，同时保留 +/- epsilon 带
yMaxA = max(E_taylor, targetError * 3) * 1.15;
ylim(axA, [-yMaxA * 0.6, yMaxA]);
xlabel(axA, '归一化波数 kh / \\pi');
ylabel(axA, '有符号绝对误差 k_{num}\\Delta - k\\Delta');
title(axA, sprintf('通带内频散误差（kh/\\pi \\in [0, %.2f]）', xMaxA), 'FontSize', 12);
legend(axA, 'Location', 'northeast', 'FontSize', 10, 'Box', 'on');
grid(axA, 'on'); box(axA, 'on');
set(axA, 'LineWidth', 0.9, 'FontSize', 11);

% --- 126 倍数值标注框 ---
annotationStr = sprintf(['\\bf 通带内 (kh \\leq %.2f\\pi) 最大绝对误差:\\rm\n', ...
    'Taylor  %d阶:  max|k_{num}\\Delta - k\\Delta| = %.4e\n', ...
    '优化    %d阶:  max|k_{num}\\Delta - k\\Delta| = %.4e\n', ...
    '\\bf 改善倍数:  %.4e / %.4e = %.1f\\times\\rm'], ...
    xCut, 2*order, E_taylor, 2*order, E_opt, E_taylor, E_opt, improvement);
annotation(figA, 'textbox', [0.42, 0.28, 0.38, 0.18], ...
    'String', annotationStr, ...
    'BackgroundColor', [1 1 0.95], 'EdgeColor', [0.2 0.2 0.2], ...
    'LineWidth', 1.0, 'FontSize', 10, 'Interpreter', 'tex', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FitBoxToText', 'off');

% 保存
outDir = fullfile(projectDir, 'results_standard', 'figures');
ewm_ensure_dir(outDir);
outFileA = fullfile(outDir, 'dispersion_passband_equiripple.png');
ewm_save_figure(figA, outFileA);
fprintf('图A 已保存: %s\n', outFileA);
close(figA);

% ======================================================================
% 4. 图B：全波数域绝对误差（对数轴）
% ======================================================================
fprintf('\n=== 生成图B：全波数域绝对误差图 ===\n');

absErrTaylor    = abs(errTaylor);
absErrOptimized = abs(errOptimized);
% 对数轴下限裁剪，避免 log(0)
floorVal = max(targetError * 1e-3, 1e-10);
absErrTaylorPlot    = max(absErrTaylor,    floorVal);
absErrOptimizedPlot = max(absErrOptimized, floorVal);

figB = figure('Color', 'w', 'Position', [100, 100, 1000, 700]);
tloB = tiledlayout(figB, 1, 1, 'Padding', 'compact');
tloB.Title.String = sprintf(['全波数域绝对误差对比 (%d阶交错网格, ' ...
    '\\epsilon = 10^{-4}, 改善 %.1f\\times)'], 2*order, improvement);
tloB.Title.FontWeight = 'bold';
tloB.Title.FontSize = 13;
tloB.Title.Interpreter = 'tex';

axB = nexttile(1); hold(axB, 'on');

% Taylor 曲线
semilogy(axB, xNorm, absErrTaylorPlot, '-', 'Color', [0.88, 0.25, 0.05], ...
    'LineWidth', 1.6, 'DisplayName', sprintf('常规 %d 阶 (Taylor)', 2*order));

% 优化曲线
semilogy(axB, xNorm, absErrOptimizedPlot, '-', 'Color', [0.00, 0.36, 0.90], ...
    'LineWidth', 2.2, 'DisplayName', sprintf('优化 %d 阶 (SA minimax)', 2*order));

% epsilon 参考线
yline(axB, targetError, 'k--', 'LineWidth', 1.0, ...
    'DisplayName', '\\epsilon = 10^{-4}');

% 优化带宽截止线
xline(axB, xCut, '--', 'Color', [0.3, 0.3, 0.3], 'LineWidth', 1.2, ...
    'DisplayName', sprintf('优化带宽 k_c = %.2f\\pi', xCut));

% --- 区域标注（通带 / 过渡区 / 高波数发散区） ---
% 通带区
txtPassband = text(axB, xCut/2, targetError * 0.3, ...
    sprintf('通带\n(SA 等纹波压制\n改善 %.0f\\times)', improvement), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.0, 0.4, 0.7], ...
    'BackgroundColor', [0.95, 0.97, 1.0], 'EdgeColor', [0.0, 0.4, 0.7], ...
    'LineWidth', 0.8, 'Margin', 4);

% 过渡区
txtTrans = text(axB, (xCut + 0.78)/2, targetError * 5, ...
    '过渡区', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.5, 0.3, 0.0], ...
    'BackgroundColor', [1.0, 0.97, 0.90], 'EdgeColor', [0.6, 0.4, 0.0], ...
    'LineWidth', 0.8, 'Margin', 4);

% 高波数发散区
txtHigh = text(axB, 0.90, targetError * 2000, ...
    '高波数发散区', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.7, 0.1, 0.1], ...
    'BackgroundColor', [1.0, 0.93, 0.93], 'EdgeColor', [0.7, 0.1, 0.1], ...
    'LineWidth', 0.8, 'Margin', 4);

% --- 坐标轴 ---
xlim(axB, [0, 1.0]);
ylim(axB, [floorVal, max([absErrTaylorPlot; absErrOptimizedPlot]) * 3]);
xlabel(axB, '归一化波数 kh / \\pi');
ylabel(axB, '绝对误差 |k_{num}\\Delta - k\\Delta|');
title(axB, '全 Nyquist 范围频散绝对误差（对数轴）', 'FontSize', 12);
legend(axB, 'Location', 'southeast', 'FontSize', 10, 'Box', 'on');
grid(axB, 'on'); box(axB, 'on');
set(axB, 'LineWidth', 0.9, 'FontSize', 11);

% --- 126 倍数值标注框 ---
annotationB = sprintf(['\\bf 通带内 (kh \\leq %.2f\\pi) 最大绝对误差:\\rm\n', ...
    'Taylor  %d阶:  max|err| = %.4e\n', ...
    '优化    %d阶:  max|err| = %.4e\n', ...
    '\\bf 改善倍数:  %.4e / %.4e = %.1f\\times\\rm'], ...
    xCut, 2*order, E_taylor, 2*order, E_opt, E_taylor, E_opt, improvement);
annotation(figB, 'textbox', [0.52, 0.20, 0.36, 0.16], ...
    'String', annotationB, ...
    'BackgroundColor', [1 1 0.95], 'EdgeColor', [0.2 0.2 0.2], ...
    'LineWidth', 1.0, 'FontSize', 10, 'Interpreter', 'tex', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FitBoxToText', 'off');

% 保存
outFileB = fullfile(outDir, 'dispersion_fullband_logerror.png');
ewm_save_figure(figB, outFileB);
fprintf('图B 已保存: %s\n', outFileB);
close(figB);

% ======================================================================
fprintf('\n=== 全部完成 ===\n');
fprintf('图A (通带等纹波): %s\n', outFileA);
fprintf('图B (全波数对数): %s\n', outFileB);
fprintf('\n关键数据:\n');
fprintf('  E_taylor = %.15g\n', E_taylor);
fprintf('  E_opt    = %.15g\n', E_opt);
fprintf('  改善倍数 = %.3f\n', improvement);


% ====================== 局部函数 ======================

function err = staggered_signed_error(coeff, kh)
% 交错网格一阶导数的有符号绝对误差 (SA 目标函数形式):
%   k_num*Delta = 2 * sum_{n=0..N-1} c_{n+1} * sin((n+0.5)*kh)
%   err         = k_num*Delta - kh
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
    error('未找到系数文件: %s', coeffFile);
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
    error('coefficients.txt 中未找到 Taylor 或优化系数。');
end
end

function s = ternary(cond, t, f)
if cond, s = t; else, s = f; end
end
