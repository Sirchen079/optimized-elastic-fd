function ewm_plot_comparison(model, leftResult, rightResult, leftLabel, rightLabel, outFile, fixedLim)
%EWM_PLOT_COMPARISON 保存并排快照对比图。

if nargin < 7
    fixedLim = [];
end

ewm_apply_chinese_style();
left = leftResult.snapshots.vz;
right = rightResult.snapshots.vz;
nSnap = min(size(left, 3), size(right, 3));
times = leftResult.snapshots.time(1:nSnap);

if isempty(fixedLim)
    commonLim = robust_wave_limit(left(:, :, 1:nSnap), right(:, :, 1:nSnap));
else
    commonLim = fixedLim;
end

% 裁掉左右空白：只显示存在显著波场的水平区段。
xKm = model.x / 1000;
zKm = model.z / 1000;
xl = ewm_wave_xlim(xKm, left(:, :, 1:nSnap), right(:, :, 1:nSnap));

% 按裁剪后波场的宽高比设置画布尺寸，使 axis image 既不留左右空白、
% 也不留上下空白；同时整体放大画布以容纳更大的坐标轴与标题字号。
% 注：标题与 xlabel 会占去行高，真实绘图区高度比 rowHeight 小，
% 因此用 effHeight 估算面板实际宽度，避免画布过宽、两列之间留大片空白。
panelAspect = (xl(2) - xl(1)) / (zKm(end) - zKm(1));   % 面板宽/高
rowHeight = 440;                                        % 每行波场的绘制高度（px）
effHeight = rowHeight - 84;                             % 扣除标题与 xlabel 后的实际绘图高度
panelWidth = effHeight * panelAspect;
figWidth = round(2 * panelWidth + 2 * 132);             % 每列仅预留色标+ylabel边距
figHeight = round(nSnap * (rowHeight + 92) + 50);       % 每行预留标题+xlabel边距
figWidth = min(max(figWidth, 760), 2400);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, figWidth, figHeight]);
tiledlayout(fig, nSnap, 2, 'Padding', 'tight', 'TileSpacing', 'compact');

for k = 1:nSnap
    lim = commonLim;

    nexttile;
    imagesc(xKm, zKm, left(:, :, k));
    axis image;
    xlim(xl);
    set(gca, 'YDir', 'reverse');
    colormap(gca, ewm_wavefield_colormap());
    caxis([-lim, lim]);
    colorbar;
    title(sprintf('%s，t = %.3f s', leftLabel, times(k)));
    xlabel('水平距离 (km)');
    ylabel('深度 (km)');

    nexttile;
    imagesc(xKm, zKm, right(:, :, k));
    axis image;
    xlim(xl);
    set(gca, 'YDir', 'reverse');
    colormap(gca, ewm_wavefield_colormap());
    caxis([-lim, lim]);
    colorbar;
    title(sprintf('%s，t = %.3f s', rightLabel, times(k)));
    xlabel('水平距离 (km)');
    ylabel('深度 (km)');
end

ewm_save_figure(fig, outFile);
close(fig);
end

function lim = robust_wave_limit(leftField, rightField)
% 使用截断的显示范围，使弱波前在图中仍然可见。
values = abs([leftField(:); rightField(:)]);
values = values(isfinite(values) & values > 0);
if isempty(values)
    lim = 1;
    return;
end

maxValue = max(values);
values = sort(values);
idx = max(1, min(numel(values), round(0.985 * numel(values))));
percentileValue = values(idx);
lim = min(maxValue, max(percentileValue, 0.22 * maxValue));
if lim == 0
    lim = 1;
end
end
