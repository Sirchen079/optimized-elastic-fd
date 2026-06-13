function ewm_plot_energy(leftResult, rightResult, leftLabel, rightLabel, outFile)
%EWM_PLOT_ENERGY 保存总能量与边界能量诊断图。

ewm_apply_chinese_style();

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1180, 460]);
tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

leftEnergy = leftResult.energy.totalVelocity;
rightEnergy = rightResult.energy.totalVelocity;
scale = max([leftEnergy(:); rightEnergy(:); eps]);

nexttile;
plot(leftResult.energy.time, leftEnergy / scale, 'b-', 'LineWidth', 1.8);
hold on;
plot(rightResult.energy.time, rightEnergy / scale, 'r-', 'LineWidth', 1.8);
style_axes();
xlabel('时间 (s)');
ylabel('归一化速度平方和');
title('物理模型区域内剩余波场能量');
legend({leftLabel, rightLabel}, 'Location', 'northeast');
ylim([0, 1.05]);
finalRatio = rightEnergy(end) / (leftEnergy(end) + eps);
text(0.03, 0.10, sprintf('末时刻 PML / 无吸收 = %.3f', finalRatio), ...
    'Units', 'normalized', ...
    'BackgroundColor', 'w', ...
    'EdgeColor', [0.72, 0.72, 0.72], ...
    'Margin', 6);

nexttile;
plot(leftResult.energy.time, leftResult.energy.boundaryRatio, 'b-', 'LineWidth', 1.8);
hold on;
plot(rightResult.energy.time, rightResult.energy.boundaryRatio, 'r-', 'LineWidth', 1.8);
style_axes();
xlabel('时间 (s)');
ylabel('边缘区域能量占比');
title('物理模型边缘附近波场能量占比');
legend({leftLabel, rightLabel}, 'Location', 'northwest');

ewm_save_figure(fig, outFile);
close(fig);
end

function style_axes()
grid on;
box on;
set(gca, 'LineWidth', 0.9, 'FontSize', 11, 'GridAlpha', 0.18);
end
