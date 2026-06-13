function regen_const_velocity_trace_hires()
%REGEN_CONST_VELOCITY_TRACE_HIRES 高 DPI 重画 const_velocity_trace_compare.png。
%
% 直接从 results_const_velocity 下已保存的 .mat 中取出 trace，
% 不重跑常速度模拟。

projectDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectDir, 'src'));
ewm_apply_chinese_style();

outDir = fullfile(projectDir, 'results_const_velocity');
figDir = fullfile(outDir, 'figures');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

taylorFile = fullfile(outDir, 'const_velocity_taylor.mat');
optFile    = fullfile(outDir, 'const_velocity_optimized.mat');
if ~isfile(taylorFile) || ~isfile(optFile)
    error('regen:MissingResults', '缺少 const_velocity 结果文件。请先跑 run_const_velocity_coeff_comparison.m。');
end
T = load(taylorFile);
O = load(optFile);
resultTaylor = T.resultTaylor;
resultOpt    = O.resultOpt;

% 接收点坐标从 summary.mat 中读，缺失则回退默认。
summaryFile = fullfile(outDir, 'summary.mat');
recZkm = 3.0; recXkm = 5.0;
if isfile(summaryFile)
    S = load(summaryFile);
    if isfield(S, 'results') && isfield(S.results, 'receiver')
        recZkm = S.results.receiver.zKm;
        recXkm = S.results.receiver.xKm;
    end
end

outFile = fullfile(figDir, 'const_velocity_trace_compare.png');
plot_trace_compare(resultTaylor, resultOpt, ...
    'Taylor 系数', '基于最大范数目标函数的优化系数', ...
    [recZkm, recXkm], outFile);
fprintf('已重新生成：%s\n', outFile);
end

function plot_trace_compare(leftResult, rightResult, ...
    leftLabel, rightLabel, recZX, outFile)
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
