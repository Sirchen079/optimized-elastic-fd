function zoomInfo = ewm_plot_wavefield_zoom(model, leftResult, rightResult, leftLabel, rightLabel, outFile, infoFile)
%EWM_PLOT_WAVEFIELD_ZOOM 保存最大差异区域的局部放大图。

if nargin < 7
    infoFile = '';
end

ewm_apply_chinese_style();

left = leftResult.snapshots.vz;
right = rightResult.snapshots.vz;
nSnap = min(size(left, 3), size(right, 3));
snapIndex = nSnap;
time = leftResult.snapshots.time(snapIndex);

leftField = left(:, :, snapIndex);
rightField = right(:, :, snapIndex);
diffField = rightField - leftField;

[~, linearIndex] = max(abs(diffField(:)));
[centerZ, centerX] = ind2sub(size(diffField), linearIndex);
if isempty(centerZ) || centerZ < 1
    centerZ = round(model.nz / 2);
    centerX = round(model.nx / 2);
end

halfZ = max(4, round(0.18 * model.nz));
halfX = max(4, round(0.18 * model.nx));
zRange = max(1, centerZ - halfZ):min(model.nz, centerZ + halfZ);
xRange = max(1, centerX - halfX):min(model.nx, centerX + halfX);

waveLim = robust_limit(abs([leftField(zRange, xRange); rightField(zRange, xRange)]), 0.99, 0.15);
diffField_zoom = diffField(zRange, xRange);
diffLim = robust_limit(abs(diffField_zoom), 0.99, 0.08);
diffPeak = max(abs(diffField_zoom(:)));

zoomWidthKm = (model.x(xRange(end)) - model.x(xRange(1))) / 1000;
zoomHeightKm = (model.z(zRange(end)) - model.z(zRange(1))) / 1000;
aspectWH = max(zoomWidthKm / max(zoomHeightKm, eps), 1);

if aspectWH >= 4
    figLayout = 'stack';
    figWidth = 1280;
    figHeight = 760;
else
    figLayout = 'sideBySide';
    figWidth = 1560;
    figHeight = 520;
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, figWidth, figHeight]);

annotation(fig, 'textbox', [0.06, 0.945, 0.88, 0.045], ...
    'String', sprintf('实验 3 波场差异局部放大区域：t = %.3f s（深度 %.2f–%.2f km，水平 %.2f–%.2f km）', ...
        time, model.z(zRange(1))/1000, model.z(zRange(end))/1000, ...
        model.x(xRange(1))/1000, model.x(xRange(end))/1000), ...
    'HorizontalAlignment', 'center', 'EdgeColor', 'none', ...
    'FontWeight', 'bold', 'FontSize', 13);

switch figLayout
    case 'stack'
        tl = tiledlayout(fig, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    case 'sideBySide'
        tl = tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
end
tl.OuterPosition = [0, 0, 1, 0.93];

nexttile;
plot_field(model, xRange, zRange, leftField(zRange, xRange), waveLim, ...
    sprintf('%s局部放大', leftLabel), 'm/s');

nexttile;
plot_field(model, xRange, zRange, rightField(zRange, xRange), waveLim, ...
    sprintf('%s局部放大', rightLabel), 'm/s');

nexttile;
plot_field(model, xRange, zRange, diffField_zoom, diffLim, ...
    sprintf('差值场局部放大（峰值 |\\Delta| = %.3g）', diffPeak), 'm/s');

ewm_save_figure(fig, outFile);
close(fig);

zoomInfo = struct();
zoomInfo.time = time;
zoomInfo.centerIndex = [centerZ, centerX];
zoomInfo.zIndexRange = [zRange(1), zRange(end)];
zoomInfo.xIndexRange = [xRange(1), xRange(end)];
zoomInfo.zKmRange = [model.z(zRange(1)), model.z(zRange(end))] / 1000;
zoomInfo.xKmRange = [model.x(xRange(1)), model.x(xRange(end))] / 1000;

if ~isempty(infoFile)
    fid = fopen(infoFile, 'w');
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '局部放大时刻_s = %.15g\n', zoomInfo.time);
    fprintf(fid, '中心z索引 = %d\n', zoomInfo.centerIndex(1));
    fprintf(fid, '中心x索引 = %d\n', zoomInfo.centerIndex(2));
    fprintf(fid, 'z索引下限 = %d\n', zoomInfo.zIndexRange(1));
    fprintf(fid, 'z索引上限 = %d\n', zoomInfo.zIndexRange(2));
    fprintf(fid, 'x索引下限 = %d\n', zoomInfo.xIndexRange(1));
    fprintf(fid, 'x索引上限 = %d\n', zoomInfo.xIndexRange(2));
    fprintf(fid, '深度下限_km = %.15g\n', zoomInfo.zKmRange(1));
    fprintf(fid, '深度上限_km = %.15g\n', zoomInfo.zKmRange(2));
    fprintf(fid, '水平距离下限_km = %.15g\n', zoomInfo.xKmRange(1));
    fprintf(fid, '水平距离上限_km = %.15g\n', zoomInfo.xKmRange(2));
end
end

function plot_field(model, xRange, zRange, field, limValue, panelTitle, cbUnit)
if nargin < 7
    cbUnit = '';
end
ax = gca;
imagesc(ax, model.x(xRange) / 1000, model.z(zRange) / 1000, field);
axis(ax, 'tight');
set(ax, 'YDir', 'reverse', ...
    'DataAspectRatioMode', 'auto', 'PlotBoxAspectRatioMode', 'auto', ...
    'Layer', 'top', 'TickDir', 'out', 'LineWidth', 0.8);
colormap(ax, ewm_wavefield_colormap());
caxis(ax, [-limValue, limValue]);
cb = colorbar(ax);
cb.Ticks = linspace(-limValue, limValue, 5);
cb.TickLabels = arrayfun(@(v) format_tick(v), cb.Ticks, 'UniformOutput', false);
if ~isempty(cbUnit)
    cb.Label.String = cbUnit;
end
cb.TickLength = 0.02;
cb.LineWidth = 0.8;
xlabel(ax, '水平距离 (km)');
ylabel(ax, '深度 (km)');
title(ax, panelTitle);
end

function s = format_tick(v)
if v == 0
    s = '0';
elseif abs(v) >= 0.1 && abs(v) < 1000
    s = sprintf('%.2f', v);
else
    s = sprintf('%.2g', v);
end
end

function lim = robust_limit(values, percentile, floorFraction)
values = values(isfinite(values) & values > 0);
if isempty(values)
    lim = 1;
    return;
end
maxValue = max(values);
values = sort(values);
idx = max(1, min(numel(values), round(percentile * numel(values))));
lim = min(maxValue, max(values(idx), floorFraction * maxValue));
if lim == 0
    lim = 1;
end
end
