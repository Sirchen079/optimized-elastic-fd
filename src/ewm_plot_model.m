function ewm_plot_model(model, cfg)
%EWM_PLOT_MODEL 保存 Marmousi 模型 Vp、Vs、密度面板图。

ewm_apply_chinese_style();
% 数据宽高比 ≈ 17:3.5 ≈ 4.86，每个子图自然又宽又扁。固定画布会让
% axis image 把扁面板在高格子里居中，从而在面板上下留出大片空白。
% 这里按数据宽高比反推面板实际高度，使每行格子高度≈面板高度，去掉竖向空白。
xKm = model.x / 1000;
zKm = model.z / 1000;
dataAspect = (xKm(end) - xKm(1)) / (zKm(end) - zKm(1));   % 宽/高
figWidth = 1360;
panelWidth = figWidth - 178;          % 扣除左侧 ylabel 与右侧色标边距
panelHeight = panelWidth / dataAspect;
rowExtra = 96;                        % 每行标题 + xlabel + 间距所需高度
figHeight = round(3 * (panelHeight + rowExtra) + 24);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, figWidth, figHeight]);
tl = tiledlayout(fig, 3, 1, 'Padding', 'tight', 'TileSpacing', 'tight');

panelLabels = {'(a)', '(b)', '(c)'};
panelTitles = {'P 波速度 Vp (km/s)', 'S 波速度 Vs (km/s)', '密度 (g/cm^3)'};
panelData = {model.vp / 1000, model.vs / 1000, model.rho / 1000};

for k = 1:3
    ax = nexttile(tl);
    imagesc(ax, model.x / 1000, model.z / 1000, panelData{k});
    axis(ax, 'image');
    set(ax, 'YDir', 'reverse');
    title(ax, panelTitles{k});
    xlabel(ax, '水平距离 (km)');
    ylabel(ax, '深度 (km)');
    colorbar(ax);
    text(ax, 0.012, 0.94, panelLabels{k}, 'Units', 'normalized', ...
        'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k', ...
        'BackgroundColor', [1 1 1 0.75], 'Margin', 2, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end

ewm_save_figure(fig, fullfile(cfg.output.dir, 'figures', 'marmousi_model.png'));
close(fig);
write_model_stats(fullfile(cfg.output.dir, 'marmousi_model_stats.csv'), model);

% --- 专用 Vs 单图（更大尺寸，清晰展示修改后的 S 波速度分布）---
vsOffset = 0;
if isfield(cfg, 'model') && isfield(cfg.model, 'vsOffset')
    vsOffset = cfg.model.vsOffset;
end
plot_vs_standalone(model, vsOffset, fullfile(cfg.output.dir, 'figures', 'vs_model.png'));
end

function plot_vs_standalone(model, vsOffset, outFile)
ewm_apply_chinese_style();
% 同样按数据宽高比反推画布高度，避免扁平面板上下留大量空白。
xKm = model.x / 1000;
zKm = model.z / 1000;
dataAspect = (xKm(end) - xKm(1)) / (zKm(end) - zKm(1));
figWidth = 1360;
panelHeight = (figWidth - 178) / dataAspect;
figHeight = round(panelHeight + 132);   % 预留标题与 xlabel 高度
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, figWidth, figHeight]);

ax = axes(fig);
vsKm = model.vs / 1000;
imagesc(ax, model.x / 1000, model.z / 1000, vsKm);
axis(ax, 'image');
set(ax, 'YDir', 'reverse');
colormap(ax, parula);
cb = colorbar(ax);
cb.Label.String = 'Vs (km/s)';
cb.Label.FontSize = 16;
xlabel(ax, '水平距离 (km)');
ylabel(ax, '深度 (km)');

vsMin = min(vsKm(:));
vsMax = max(vsKm(:));
if vsOffset ~= 0
    titleStr = sprintf('S 波速度 Vs（在原始模型基础上 +%g m/s）    实际范围：%.3f – %.3f km/s', ...
        vsOffset, vsMin, vsMax);
else
    titleStr = sprintf('S 波速度 Vs    范围：%.3f – %.3f km/s', vsMin, vsMax);
end
title(ax, titleStr);

ewm_save_figure(fig, outFile);
close(fig);
end

function write_model_stats(outFile, model)
fid = fopen(outFile, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '参数,网格点数_z,网格点数_x,dz_m,dx_m,最小值,最大值,平均值,单位\n');
write_row(fid, 'Vp', model, model.vp, 'm/s');
write_row(fid, 'Vs', model, model.vs, 'm/s');
write_row(fid, 'rho', model, model.rho, 'kg/m^3');
end

function write_row(fid, name, model, values, unitText)
fprintf(fid, '%s,%d,%d,%.15g,%.15g,%.15g,%.15g,%.15g,%s\n', ...
    name, model.nz, model.nx, model.dz, model.dx, ...
    min(values(:)), max(values(:)), mean(values(:)), unitText);
end
