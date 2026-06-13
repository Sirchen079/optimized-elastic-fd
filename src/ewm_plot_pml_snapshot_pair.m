function ewm_plot_pml_snapshot_pair(model, noAbsorb, pml, outFile)
%EWM_PLOT_PML_SNAPSHOT_PAIR PML 对比的单时刻快照对。

ewm_apply_chinese_style();

left = noAbsorb.snapshots.vz(:, :, end);
right = pml.snapshots.vz(:, :, end);
time = noAbsorb.snapshots.time(end);

lim = 2;

% 裁掉左右空白：只显示存在显著波场的水平区段，并按其宽高比设置画布。
xKm = model.x / 1000;
zKm = model.z / 1000;
xl = ewm_wave_xlim(xKm, left, right);
panelAspect = (xl(2) - xl(1)) / (zKm(end) - zKm(1));
rowHeight = 380;
panelWidth = rowHeight * panelAspect;
figWidth = round(2 * panelWidth + 2 * 175 + 90);
figHeight = round(rowHeight + 130);
figWidth = min(max(figWidth, 980), 2400);

% 本图为单行布局、整体偏宽，放进论文按页宽缩放后字号会显得偏小。
% 这里按图宽相对参考宽度等比放大字号，使其视觉大小与多行快照图一致。
fontScale = max(1, figWidth / 1180);
axFont = round(18 * fontScale);
cbFont = round(16 * fontScale);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, figWidth, figHeight]);
tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
imagesc(xKm, zKm, left);
axis image;
xlim(xl);
set(gca, 'YDir', 'reverse', 'FontSize', axFont);
colormap(gca, ewm_wavefield_colormap());
caxis([-lim, lim]);
cb1 = colorbar;
cb1.FontSize = cbFont;
title(sprintf('无吸收边界，t = %.3f s', time));
xlabel('水平距离 (km)');
ylabel('深度 (km)');

nexttile;
imagesc(xKm, zKm, right);
axis image;
xlim(xl);
set(gca, 'YDir', 'reverse', 'FontSize', axFont);
colormap(gca, ewm_wavefield_colormap());
caxis([-lim, lim]);
cb2 = colorbar;
cb2.FontSize = cbFont;
title(sprintf('PML 吸收边界，t = %.3f s', time));
xlabel('水平距离 (km)');
ylabel('深度 (km)');

ewm_save_figure(fig, outFile);
close(fig);
end
