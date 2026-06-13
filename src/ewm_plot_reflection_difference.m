function ewm_plot_reflection_difference(model, noAbsorb, pml, outFile)
%EWM_PLOT_REFLECTION_DIFFERENCE 展示边界反射残差。

ewm_apply_chinese_style();

left = noAbsorb.snapshots.vz;
right = pml.snapshots.vz;
nSnap = min(size(left, 3), size(right, 3));
times = noAbsorb.snapshots.time(1:nSnap);

figHeight = max(420, 260 * nSnap);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1460, figHeight]);
tiledlayout(fig, nSnap, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for k = 1:nSnap
    diffField = left(:, :, k) - right(:, :, k);
    waveLim = max(abs([left(:, :, k); right(:, :, k)]), [], 'all');
    diffLim = max(abs(diffField), [], 'all');
    if waveLim == 0
        waveLim = 1;
    end
    if diffLim == 0
        diffLim = 1;
    end

    nexttile;
    imagesc(model.x / 1000, model.z / 1000, left(:, :, k));
    axis image;
    set(gca, 'YDir', 'reverse');
    colormap(gca, ewm_wavefield_colormap());
    caxis([-waveLim, waveLim]);
    colorbar;
    title(sprintf('无吸收边界，t = %.3f s', times(k)));
    xlabel('水平距离 (km)');
    ylabel('深度 (km)');

    nexttile;
    imagesc(model.x / 1000, model.z / 1000, right(:, :, k));
    axis image;
    set(gca, 'YDir', 'reverse');
    colormap(gca, ewm_wavefield_colormap());
    caxis([-waveLim, waveLim]);
    colorbar;
    title(sprintf('PML 吸收边界，t = %.3f s', times(k)));
    xlabel('水平距离 (km)');
    ylabel('深度 (km)');

    nexttile;
    imagesc(model.x / 1000, model.z / 1000, diffField);
    axis image;
    set(gca, 'YDir', 'reverse');
    colormap(gca, ewm_wavefield_colormap());
    caxis([-diffLim, diffLim]);
    colorbar;
    title(sprintf('差异场：无吸收 - PML，t = %.3f s', times(k)));
    xlabel('水平距离 (km)');
    ylabel('深度 (km)');
end

ewm_save_figure(fig, outFile);
close(fig);
end
