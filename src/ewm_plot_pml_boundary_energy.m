function ewm_plot_pml_boundary_energy(noAbsorb, pml, outFile, csvFile)
%EWM_PLOT_PML_BOUNDARY_ENERGY 绘制边界区域速度能量随时间变化曲线。

ewm_apply_chinese_style();

n = min(numel(noAbsorb.energy.time), numel(pml.energy.time));
time = noAbsorb.energy.time(1:n);
noBoundary = boundary_energy(noAbsorb, n);
pmlBoundary = boundary_energy(pml, n);
scale = max([noBoundary(:); pmlBoundary(:); eps]);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1180, 520]);

plot(time, noBoundary / scale, '-', 'Color', [0.05, 0.32, 0.62], 'LineWidth', 1.9);
hold on;
plot(time, pmlBoundary / scale, '-', 'Color', [0.82, 0.28, 0.12], 'LineWidth', 1.9);
grid on;
box on;
xlabel('时间 (s)');
ylabel('归一化边界带速度能量');
title('PML 边界带能量随时间变化');
legend({'无吸收边界', 'PML 吸收边界'}, 'Location', 'northwest');
set(gca, 'LineWidth', 0.9, 'FontSize', 11, 'GridAlpha', 0.18);

ewm_save_figure(fig, outFile);
close(fig);

fid = fopen(csvFile, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '时间_s,无吸收边界带能量,PML边界带能量,PML相对无吸收比值\n');
for k = 1:n
    fprintf(fid, '%.15g,%.15g,%.15g,%.15g\n', ...
        time(k), noBoundary(k), pmlBoundary(k), pmlBoundary(k) / (noBoundary(k) + eps));
end
end

function values = boundary_energy(result, n)
if isfield(result.energy, 'boundaryVelocity')
    values = result.energy.boundaryVelocity(1:n);
else
    values = result.energy.totalVelocity(1:n) .* result.energy.boundaryRatio(1:n);
end
values = values(:);
end
