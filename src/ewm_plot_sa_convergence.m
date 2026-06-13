function ewm_plot_sa_convergence(optInfo, targetError, outFile)
%EWM_PLOT_SA_CONVERGENCE 绘制模拟退火算法的收敛曲线。
%
% 输入：
%   optInfo      由 ewm_optimize_minimax_coeffs 返回的 info 结构，
%                必须包含 info.trace 字段。
%   targetError  目标误差阈值（用于在图上标出参考线）。
%   outFile      输出 PNG 路径。

if ~isfield(optInfo, 'trace') || isempty(optInfo.trace)
    warning('ewm:NoTrace', '优化信息中未找到 trace 字段，跳过收敛曲线绘制。');
    return;
end

trace = optInfo.trace;
step = trace.step(:);
restart = trace.restart(:);
temperature = trace.temperature(:);
bestObj = trace.bestObjective(:);
currentObj = trace.currentObjective(:);

isPolish = temperature < 0;
saMask = ~isPolish;

xAxis = (1:numel(step)).';

ewm_apply_chinese_style();
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1320, 760]);
layout = tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

stdColor = [0.05, 0.32, 0.62];
optColor = [0.82, 0.28, 0.12];
polishColor = [0.20, 0.55, 0.25];
restartColor = [0.45, 0.45, 0.45];

nexttile;
plot(xAxis(saMask), temperature(saMask), '-', 'Color', stdColor, 'LineWidth', 1.6);
hold on;
restartIdx = find_restart_boundaries(restart);
for k = 1:numel(restartIdx)
    xline(restartIdx(k) - 0.5, ':', 'Color', restartColor, 'LineWidth', 0.8, 'Alpha', 0.6);
end
set(gca, 'YScale', 'log');
grid on;
box on;
xlabel('SA 温度迭代序号');
ylabel('温度 T');
title('（a）模拟退火温度衰减曲线（虚线表示 restart 起点）');
set(gca, 'LineWidth', 0.9);

nexttile;
plot(xAxis(saMask), currentObj(saMask), '-', 'Color', [0.55 0.65 0.85], 'LineWidth', 1.0);
hold on;
plot(xAxis, bestObj, '-', 'Color', optColor, 'LineWidth', 2.2);
if any(isPolish)
    plot(xAxis(isPolish), bestObj(isPolish), 'p', 'Color', polishColor, ...
        'MarkerFaceColor', polishColor, 'MarkerSize', 10, 'LineWidth', 0.8);
end
yline(targetError, 'k--', 'LineWidth', 1.2, ...
    'Label', sprintf('目标误差 %.0e', targetError), ...
    'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');
for k = 1:numel(restartIdx)
    xline(restartIdx(k) - 0.5, ':', 'Color', restartColor, 'LineWidth', 0.8, 'Alpha', 0.6);
end
set(gca, 'YScale', 'log');
grid on;
box on;
xlabel('SA 温度迭代序号');
ylabel('最大绝对频散误差 |k_{num}\Delta - k\Delta|');
title('（b）目标函数收敛过程');
legendEntries = {'当前解目标值', '历史最优'};
if any(isPolish)
    legendEntries{end + 1} = '局部精修';
end
legendEntries{end + 1} = '目标阈值';
legend(legendEntries, 'Location', 'northeast');
set(gca, 'LineWidth', 0.9);

sgtitle(layout, '模拟退火算法收敛过程', 'FontWeight', 'bold', 'FontSize', 22);

ewm_save_figure(fig, outFile);
close(fig);
end

function idx = find_restart_boundaries(restart)
diffs = diff(restart);
boundary = find(diffs > 0) + 1;
idx = boundary(:);
end
