% 在指定空间窗口（x: 5–11 km, z: 0–1.2 km）绘制 t = 1.250 s 时
% Taylor 系数、优化系数及差值场的局部放大图。
projectDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectDir, 'src'));
ewm_apply_chinese_style();

cfg = ewm_default_config(projectDir, 'standard');
cfg.output.figuresDir = fullfile(cfg.output.dir, 'figures');

model = ewm_load_marmousi(cfg.model);

stdData = load(fullfile(cfg.output.dir, 'exp3_staggered_pml_standard.mat'));
optData = load(fullfile(cfg.output.dir, 'exp3_staggered_pml_minimax.mat'));

targetTime = 1.250;
times = stdData.result.snapshots.time;
[~, snapIndex] = min(abs(times - targetTime));
actualTime = times(snapIndex);

xLimKm = [5, 11];
zLimKm = [0, 1.2];
xRange = find(model.x >= xLimKm(1) * 1000 & model.x <= xLimKm(2) * 1000);
zRange = find(model.z >= zLimKm(1) * 1000 & model.z <= zLimKm(2) * 1000);

leftField  = stdData.result.snapshots.vz(zRange, xRange, snapIndex);
rightField = optData.result.snapshots.vz(zRange, xRange, snapIndex);
diffField  = rightField - leftField;

diffPeak = max(abs(diffField(:)));
waveLim = 2;
diffLim = 0.05;
magnification = waveLim / diffLim;

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1600, 480]);
tl = tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot_field(model, xRange, zRange, leftField, waveLim, ...
    sprintf('Taylor 系数局部放大，t = %.3f s', actualTime));

nexttile;
plot_field(model, xRange, zRange, rightField, waveLim, ...
    sprintf('基于最大范数目标函数的优化系数局部放大，t = %.3f s', actualTime));

nexttile;
plot_field(model, xRange, zRange, diffField, diffLim, '差值场局部放大');

outFile = fullfile(cfg.output.figuresDir, 'exp3_zoom_source_region_5_11km.png');
ewm_save_figure(fig, outFile);
close(fig);

fprintf('已保存：%s\n', outFile);
fprintf('实际时间 = %.6f s（请求 %.3f s）\n', actualTime, targetTime);
fprintf('水平索引 %d:%d，深度索引 %d:%d\n', xRange(1), xRange(end), zRange(1), zRange(end));
fprintf('实际窗口：x = %.3f–%.3f km，z = %.3f–%.3f km\n', ...
    model.x(xRange(1))/1000, model.x(xRange(end))/1000, ...
    model.z(zRange(1))/1000, model.z(zRange(end))/1000);
fprintf('波场色标范围 ±%.4g m/s；差值色标范围 ±%.4g m/s；差值峰值 %.4g m/s\n', ...
    waveLim, diffLim, diffPeak);

function plot_field(model, xRange, zRange, field, limValue, panelTitle)
ax = gca;
imagesc(ax, model.x(xRange) / 1000, model.z(zRange) / 1000, field);
axis(ax, 'image');
set(ax, 'YDir', 'reverse', 'Layer', 'top');
colormap(ax, ewm_wavefield_colormap());
caxis(ax, [-limValue, limValue]);
cb = colorbar(ax);
cb.Ticks = [-limValue, 0, limValue];
cb.TickLabels = {format_tick(-limValue), '0', format_tick(limValue)};
xlabel(ax, '水平距离 (km)');
ylabel(ax, '深度 (km)');
title(ax, panelTitle);
end

function s = format_tick(v)
if v == 0
    s = '0';
elseif abs(v) >= 1
    s = sprintf('%g', v);
else
    s = sprintf('%g', v);
end
end

function v = nice_ceil(x)
if x <= 0
    v = 1;
    return;
end
e = floor(log10(x));
base = x / 10^e;
nice = [1, 1.5, 2, 2.5, 3, 4, 5, 7, 10];
i = find(nice >= base - 1e-12, 1, 'first');
v = nice(i) * 10^e;
end

function lim = robust_limit(values, percentile, floorFraction)
values = values(isfinite(values) & values > 0);
if isempty(values)
    lim = 1; return;
end
maxValue = max(values);
values = sort(values(:));
idx = max(1, min(numel(values), round(percentile * numel(values))));
lim = min(maxValue, max(values(idx), floorFraction * maxValue));
if lim <= 0
    lim = 1;
end
end
