function ewm_plot_pml_reflection_comparison(model, noAbsorb, pml, outFile)
%EWM_PLOT_PML_REFLECTION_COMPARISON PML 反射直接对比。
%
% 第1、2列使用相同色标。第3列为差值场（无吸收 - PML），
% 使用独立色标。该差值场为辅助对比，并非严格的反射波场。

ewm_apply_chinese_style();

left = noAbsorb.snapshots.vz;
right = pml.snapshots.vz;
nSnap = min(size(left, 3), size(right, 3));
times = noAbsorb.snapshots.time(1:nSnap);

figHeight = max(520, 260 * nSnap);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1560, figHeight]);
tiledlayout(fig, nSnap, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for k = 1:nSnap
    noAbsorbField = left(:, :, k);
    pmlField = right(:, :, k);
    residual = noAbsorbField - pmlField;

    waveLim = max(abs([noAbsorbField(:); pmlField(:)]));
    residualLim = max(abs(residual(:)));
    if waveLim == 0
        waveLim = 1;
    end
    if residualLim == 0
        residualLim = 1;
    end

    nexttile;
    imagesc(model.x / 1000, model.z / 1000, noAbsorbField);
    axis image;
    set(gca, 'YDir', 'reverse');
    colormap(gca, ewm_wavefield_colormap());
    caxis([-waveLim, waveLim]);
    colorbar;
    hold on;
    draw_boundary_boxes(model);
    title(sprintf('无吸收边界：含边界反射，t = %.3f s', times(k)));
    xlabel('水平距离 (km)');
    ylabel('深度 (km)');

    nexttile;
    imagesc(model.x / 1000, model.z / 1000, pmlField);
    axis image;
    set(gca, 'YDir', 'reverse');
    colormap(gca, ewm_wavefield_colormap());
    caxis([-waveLim, waveLim]);
    colorbar;
    hold on;
    draw_boundary_boxes(model);
    title(sprintf('PML 吸收边界：反射被削弱，t = %.3f s', times(k)));
    xlabel('水平距离 (km)');
    ylabel('深度 (km)');

    nexttile;
    imagesc(model.x / 1000, model.z / 1000, residual);
    axis image;
    set(gca, 'YDir', 'reverse');
    colormap(gca, ewm_wavefield_colormap());
    caxis([-residualLim, residualLim]);
    colorbar;
    hold on;
    draw_boundary_boxes(model);
    title(sprintf('差异场：无吸收 - PML，t = %.3f s', times(k)));
    xlabel('水平距离 (km)');
    ylabel('深度 (km)');
end

annotation(fig, 'textbox', [0.17, 0.965, 0.70, 0.03], ...
    'String', '判读方法：前两列使用同一色标；第三列只是无吸收与 PML 的差异场，不等于纯反射场。严格判断请看“扩大计算域波场偏离”图（exp2_boundary_reference_error.png）。红框标出最容易出现边界反射的区域。', ...
    'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'FontWeight', 'bold');

ewm_save_figure(fig, outFile);
close(fig);
end

function draw_boundary_boxes(model)
x0 = model.x(1) / 1000;
x1 = model.x(end) / 1000;
z0 = model.z(1) / 1000;
z1 = model.z(end) / 1000;
width = x1 - x0;
height = z1 - z0;
boxZ = 0.18 * height;
boxX = 0.08 * width;
rectangle('Position', [x0, z0, width, boxZ], ...
    'EdgeColor', [0.85 0.1 0.1], 'LineStyle', '--', 'LineWidth', 1.2);
rectangle('Position', [x0, z1 - boxZ, width, boxZ], ...
    'EdgeColor', [0.85 0.1 0.1], 'LineStyle', '--', 'LineWidth', 1.2);
rectangle('Position', [x0, z0, boxX, height], ...
    'EdgeColor', [0.85 0.1 0.1], 'LineStyle', '--', 'LineWidth', 1.2);
rectangle('Position', [x1 - boxX, z0, boxX, height], ...
    'EdgeColor', [0.85 0.1 0.1], 'LineStyle', '--', 'LineWidth', 1.2);
end
