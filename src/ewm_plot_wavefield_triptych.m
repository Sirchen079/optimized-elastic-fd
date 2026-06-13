function metrics = ewm_plot_wavefield_triptych(model, leftResult, rightResult, leftLabel, rightLabel, outFile, metricsFile)
%EWM_PLOT_WAVEFIELD_TRIPTYCH 保存左/右/差异波场三联对比图。

if nargin < 7
    metricsFile = '';
end

ewm_apply_chinese_style();

left = leftResult.snapshots.vz;
right = rightResult.snapshots.vz;
nSnap = min(size(left, 3), size(right, 3));
times = leftResult.snapshots.time(1:nSnap);

waveLim = robust_wave_limit(left(:, :, 1:nSnap), right(:, :, 1:nSnap));
diffLim = robust_diff_limit(left(:, :, 1:nSnap), right(:, :, 1:nSnap));

relL2 = zeros(nSnap, 1);
maxAbsDiff = zeros(nSnap, 1);

figHeight = max(460, 280 * nSnap);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1560, figHeight]);
tiledlayout(fig, nSnap, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for k = 1:nSnap
    leftField = left(:, :, k);
    rightField = right(:, :, k);
    diffField = rightField - leftField;
    relL2(k) = norm(diffField(:)) / (norm(leftField(:)) + eps);
    maxAbsDiff(k) = max(abs(diffField(:)));

    nexttile;
    plot_field(model, leftField, waveLim, sprintf('%s, t = %.3f s', leftLabel, times(k)));

    nexttile;
    plot_field(model, rightField, waveLim, sprintf('%s, t = %.3f s', rightLabel, times(k)));

    nexttile;
    plot_field(model, diffField, diffLim, sprintf('差值场，相对L2 = %.3g', relL2(k)));
end

ewm_save_figure(fig, outFile);
close(fig);

metrics = struct();
metrics.time = times(:);
metrics.relativeL2 = relL2;
metrics.maxAbsDifference = maxAbsDiff;

if ~isempty(metricsFile)
    fid = fopen(metricsFile, 'w');
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '差值定义 = 优化系数波场 - Taylor系数波场\n');
    fprintf(fid, '相对L2定义 = 范数(优化系数波场 - Taylor系数波场) / 范数(Taylor系数波场)\n\n');
    fprintf(fid, '时间_s,相对L2,最大绝对差值\n');
    for k = 1:nSnap
        fprintf(fid, '%.15g,%.15g,%.15g\n', times(k), relL2(k), maxAbsDiff(k));
    end
end
end

function plot_field(model, field, limValue, panelTitle)
imagesc(model.x / 1000, model.z / 1000, field);
axis image;
set(gca, 'YDir', 'reverse');
colormap(gca, ewm_wavefield_colormap());
caxis([-limValue, limValue]);
colorbar;
xlabel('水平距离 (km)');
ylabel('深度 (km)');
title(panelTitle);
end

function lim = robust_wave_limit(leftField, rightField)
values = abs([leftField(:); rightField(:)]);
lim = robust_limit(values, 0.985, 0.22);
end

function lim = robust_diff_limit(leftField, rightField)
values = abs(rightField(:) - leftField(:));
lim = robust_limit(values, 0.985, 0.10);
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
