function outputs = ewm_plot_sa_error_evidence(mode)
%EWM_PLOT_SA_ERROR_EVIDENCE Thesis figures for optimized vs Taylor coefficients.
%
% All quantitative error annotations in this function use only the simulated
% annealing objective:
%     |k_num * Delta - k * Delta|
% The wavefield panels are qualitative evidence and are deliberately labelled
% as wavefield amplitudes or differences, not as new error definitions.

if nargin < 1 || isempty(mode)
    mode = 'standard';
end

projectDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectDir, 'src'));
ewm_apply_chinese_style();

cfg = ewm_default_config(projectDir, mode);
cfg.output.figuresDir = fullfile(cfg.output.dir, 'figures');
ewm_ensure_dir(cfg.output.figuresDir);

model = ewm_load_marmousi(cfg.model);
[standardCoeff, optimizedCoeff] = read_saved_coefficients( ...
    fullfile(cfg.output.dir, 'coefficients.txt'), cfg.coeff.order);

khMax = cfg.coeff.khMax;
targetError = cfg.coeff.targetError;

samples = 6000;
kh = linspace(khMax / samples, pi, samples).';
signedTaylor = staggered_signed_error(standardCoeff, kh);
signedOptimized = staggered_signed_error(optimizedCoeff, kh);
absTaylor = abs(signedTaylor);
absOptimized = abs(signedOptimized);
inBand = kh <= khMax + eps;

summary = struct();
summary.order = 2 * numel(standardCoeff);
summary.khMaxOverPi = khMax / pi;
summary.targetError = targetError;
summary.taylorMaxAbsError = max(absTaylor(inBand));
summary.optimizedMaxAbsError = max(absOptimized(inBand));
summary.improvement = summary.taylorMaxAbsError / summary.optimizedMaxAbsError;

outputs = struct();
outputs.figures = struct();
outputs.figures.theoryPanel = fullfile(cfg.output.figuresDir, 'sa_error_theory_panel.png');
outputs.figures.kspaceMap = fullfile(cfg.output.figuresDir, 'sa_error_kspace_map.png');
outputs.figures.wavefieldDifference = fullfile(cfg.output.figuresDir, 'exp3_sa_wavefield_difference_clean.png');
outputs.summaryFile = fullfile(cfg.output.dir, 'sa_error_standard_evidence_summary.csv');

plot_theory_panel(kh, signedTaylor, signedOptimized, absTaylor, absOptimized, ...
    inBand, summary, targetError, outputs.figures.theoryPanel);
plot_kspace_map(kh, absTaylor, absOptimized, khMax, targetError, ...
    outputs.figures.kspaceMap);

standardData = load_required_result(fullfile(cfg.output.dir, 'exp3_staggered_pml_standard.mat'));
optimizedData = load_required_result(fullfile(cfg.output.dir, 'exp3_staggered_pml_minimax.mat'));

plot_wavefield_difference_clean(model, standardData.result, optimizedData.result, ...
    outputs.figures.wavefieldDifference);

write_summary_csv(outputs.summaryFile, summary, standardCoeff, optimizedCoeff);

fprintf('SA absolute-error evidence finished. Max error: Taylor %.6g, optimized %.6g, improvement %.3f x.\n', ...
    summary.taylorMaxAbsError, summary.optimizedMaxAbsError, summary.improvement);
end

function plot_theory_panel(kh, signedTaylor, signedOptimized, absTaylor, absOptimized, inBand, summary, targetError, outFile)
colors = standard_colors();
x = kh / pi;
xCut = summary.khMaxOverPi;

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 80, 1450, 980]);
tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
hTaylorSigned = plot(x(inBand), signedTaylor(inBand), '-', 'Color', colors.taylor, 'LineWidth', 1.55);
hold on;
hOptimizedSigned = plot(x(inBand), signedOptimized(inBand), '-', 'Color', colors.optimized, 'LineWidth', 2.0);
hTargetSigned = yline(targetError, 'k--', 'LineWidth', 1.0);
yline(-targetError, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
yline(0, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.8, 'HandleVisibility', 'off');
grid on;
box on;
xlim([0, xCut]);
yLim = max(summary.taylorMaxAbsError, targetError * 4) * 1.12;
ylim([-yLim, yLim]);
xlabel('归一化波数 kh / \pi');
ylabel('有符号 SA 绝对误差  k_{num}\Delta - k\Delta');
title('(a) 通带内有符号 SA 绝对误差');
legend([hTaylorSigned, hOptimizedSigned, hTargetSigned], ...
    {'Taylor 系数', '优化系数', '\pm 1e-4 阈值'}, 'Location', 'southwest');

nexttile;
floorVal = max(targetError * 1.0e-3, 1.0e-8);
hTaylorAbs = semilogy(x(inBand), max(absTaylor(inBand), floorVal), '-', 'Color', colors.taylor, 'LineWidth', 1.55);
hold on;
hOptimizedAbs = semilogy(x(inBand), max(absOptimized(inBand), floorVal), '-', 'Color', colors.optimized, 'LineWidth', 2.0);
hTargetAbs = yline(targetError, 'k--', 'LineWidth', 1.1);
grid on;
box on;
xlim([0, xCut]);
ylim([floorVal, max(absTaylor(inBand)) * 1.8]);
xlabel('归一化波数 kh / \pi');
ylabel('|k_{num}\Delta - k\Delta|');
title('(b) 通带内 SA 绝对误差');
legend([hTaylorAbs, hOptimizedAbs, hTargetAbs], ...
    {'Taylor 系数', '优化系数', '1e-4 阈值'}, 'Location', 'southeast');

nexttile;
hTaylorFull = semilogy(x, max(absTaylor, floorVal), '-', 'Color', colors.taylor, 'LineWidth', 1.35);
hold on;
hOptimizedFull = semilogy(x, max(absOptimized, floorVal), '-', 'Color', colors.optimized, 'LineWidth', 1.85);
hTargetFull = yline(targetError, 'k--', 'LineWidth', 1.0);
hCutFull = xline(xCut, ':', 'Color', [0.20 0.20 0.20], 'LineWidth', 1.1);
hPatch = patch([0 xCut xCut 0], [floorVal floorVal max(absTaylor)*2 max(absTaylor)*2], ...
    [0.92 0.95 1.00], 'FaceAlpha', 0.26, 'EdgeColor', 'none');
set(hPatch, 'HandleVisibility', 'off');
uistack(hPatch, 'bottom');
grid on;
box on;
xlim([0, 1]);
ylim([floorVal, max(absTaylor) * 2]);
xlabel('归一化波数 kh / \pi');
ylabel('|k_{num}\Delta - k\Delta|');
title('(c) 全 Nyquist 范围 SA 绝对误差');
legend([hTaylorFull, hOptimizedFull, hTargetFull, hCutFull], ...
    {'Taylor 系数', '优化系数', '1e-4 阈值', 'kh_{max}'}, 'Location', 'southeast');

nexttile;
barValues = [summary.taylorMaxAbsError, summary.optimizedMaxAbsError];
b = bar(1:2, barValues, 0.54);
b.FaceColor = 'flat';
b.CData = [colors.taylor; colors.optimized];
set(gca, 'YScale', 'log', 'XTick', 1:2, ...
    'XTickLabel', {'Taylor', '优化'}, 'XTickLabelRotation', 0);
yline(targetError, 'k--', 'LineWidth', 1.1);
grid on;
box on;
ylim([targetError * 0.45, summary.taylorMaxAbsError * 2.2]);
ylabel('通带最大 |k_{num}\Delta - k\Delta|');
title('(d) 最大 SA 绝对误差定量对比');
for i = 1:2
    text(i, barValues(i) * 1.18, sprintf('%.3g', barValues(i)), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 15);
end
text(1.5, summary.taylorMaxAbsError * 1.75, ...
    sprintf('改善 %.1f 倍', summary.improvement), ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 16, ...
    'Color', [0.12 0.25 0.45]);

sgtitle(sprintf('Taylor 与优化系数的 SA 绝对误差对比：%d 阶，kh \\leq %.2f\\pi，阈值 %.0e', ...
    summary.order, summary.khMaxOverPi, targetError), 'FontWeight', 'bold', 'FontSize', 22);

save_raster_figure(fig, outFile);
close(fig);
end

function plot_kspace_map(kh, absTaylor, absOptimized, khMax, targetError, outFile)
gridN = 361;
axisKh = linspace(-pi, pi, gridN);
[kx, kz] = meshgrid(axisKh, axisKh);
componentKh = max(abs(kx), abs(kz));
passMask = componentKh <= khMax + eps;

floorVal = max(targetError * 1.0e-3, 1.0e-8);
taylorMap = interp1(kh, absTaylor, componentKh, 'linear', 'extrap');
optimizedMap = interp1(kh, absOptimized, componentKh, 'linear', 'extrap');
taylorMap(~passMask) = NaN;
optimizedMap(~passMask) = NaN;

logTaylor = log10(max(taylorMap, floorVal));
logOptimized = log10(max(optimizedMap, floorVal));
cl = [log10(floorVal), log10(max(taylorMap(:), [], 'omitnan'))];

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1380, 620]);
tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
imagesc(axisKh / pi, axisKh / pi, logTaylor);
axis image;
set(gca, 'YDir', 'normal');
colormap(gca, parula(256));
caxis(cl);
hold on;
draw_passband_box(khMax / pi);
title('(a) Taylor 系数');
xlabel('k_x h / \pi');
ylabel('k_z h / \pi');
cb = colorbar;
format_log_colorbar(cb);

nexttile;
imagesc(axisKh / pi, axisKh / pi, logOptimized);
axis image;
set(gca, 'YDir', 'normal');
colormap(gca, parula(256));
caxis(cl);
hold on;
draw_passband_box(khMax / pi);
title('(b) 优化系数');
xlabel('k_x h / \pi');
ylabel('k_z h / \pi');
cb = colorbar;
format_log_colorbar(cb);

sgtitle(sprintf('SA 绝对误差在二维波数平面的通带映射：|k_{num}\\Delta - k\\Delta|，阈值 %.0e', targetError), ...
    'FontWeight', 'bold');

save_raster_figure(fig, outFile);
close(fig);
end

function plot_wavefield_difference_clean(model, standardResult, optimizedResult, outFile)
standard = standardResult.snapshots.vz;
optimized = optimizedResult.snapshots.vz;
nSnap = min(size(standard, 3), size(optimized, 3));
snapIndex = nSnap;
time = standardResult.snapshots.time(snapIndex);

waveLim = robust_limit(abs([standard(:, :, 1:nSnap); optimized(:, :, 1:nSnap)]), 0.985, 0.22);
diffLim = robust_limit(abs(optimized(:, :, 1:nSnap) - standard(:, :, 1:nSnap)), 0.985, 0.12);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 80, 1540, 520]);

stdField = standard(:, :, snapIndex);
optField = optimized(:, :, snapIndex);
diffField = optField - stdField;
[zRange, xRange] = active_wavefield_window(stdField, optField, model);

annotation(fig, 'textbox', [0.06, 0.91, 0.88, 0.06], ...
    'String', '实验 3 代表性波场定性对比（有效波场窗口）：同一模型、同一时间步长、同一 PML 条件', ...
    'HorizontalAlignment', 'center', 'EdgeColor', 'none', ...
    'FontWeight', 'bold', 'FontSize', 14);

ax1 = axes(fig, 'Position', [0.055, 0.17, 0.265, 0.68]);
plot_field_window(ax1, model, zRange, xRange, stdField(zRange, xRange), waveLim, ...
    sprintf('Taylor 波场，t = %.3f s', time));

ax2 = axes(fig, 'Position', [0.385, 0.17, 0.265, 0.68]);
plot_field_window(ax2, model, zRange, xRange, optField(zRange, xRange), waveLim, ...
    sprintf('优化系数波场，t = %.3f s', time));

ax3 = axes(fig, 'Position', [0.715, 0.17, 0.265, 0.68]);
plot_field_window(ax3, model, zRange, xRange, diffField(zRange, xRange), diffLim, ...
    '差值场：优化系数 - Taylor');

ewm_save_figure(fig, outFile);
close(fig);
end

function plot_field_window(ax, model, zRange, xRange, field, limValue, panelTitle)
imagesc(ax, model.x(xRange) / 1000, model.z(zRange) / 1000, field);
axis(ax, 'tight');
set(ax, 'DataAspectRatioMode', 'auto', 'PlotBoxAspectRatioMode', 'auto');
set(ax, 'YDir', 'reverse');
colormap(ax, ewm_wavefield_colormap());
caxis(ax, [-limValue, limValue]);
colorbar(ax);
xlabel(ax, '水平距离 (km)');
ylabel(ax, '深度 (km)');
title(ax, panelTitle);
end

function [zRange, xRange] = active_wavefield_window(fieldA, fieldB, model)
amp = max(abs(fieldA), abs(fieldB));
threshold = 0.025 * max(amp(:));
mask = amp >= threshold;
if ~any(mask(:))
    zRange = 1:model.nz;
    xRange = 1:model.nx;
    return;
end
[zIdx, xIdx] = find(mask);
padZ = max(6, round(0.10 * model.nz));
padX = max(12, round(0.08 * model.nx));
z1 = max(1, min(zIdx) - padZ);
z2 = min(model.nz, max(zIdx) + padZ);
x1 = max(1, min(xIdx) - padX);
x2 = min(model.nx, max(xIdx) + padX);
zRange = z1:z2;
xRange = x1:x2;
end

function draw_passband_box(limit)
plot([-limit limit limit -limit -limit], [-limit -limit limit limit -limit], ...
    'k-', 'LineWidth', 1.15);
end

function format_log_colorbar(cb)
ticks = -8:1:-1;
cb.Ticks = ticks;
labels = cell(size(ticks));
for i = 1:numel(ticks)
    labels{i} = sprintf('10^{%d}', ticks(i));
end
cb.TickLabels = labels;
end

function save_raster_figure(fig, outFile)
[outDir, ~, ~] = fileparts(outFile);
ewm_ensure_dir(outDir);
warningState = warning;
cleanupWarning = onCleanup(@() warning(warningState));
warning('off', 'all');
try
    fig.Renderer = 'opengl';
catch
end
drawnow;
try
    exportgraphics(fig, outFile, 'Resolution', 900);
catch
    print(fig, outFile, '-dpng', '-r900', '-opengl');
end
end

function lim = robust_limit(values, percentile, floorFraction)
values = values(isfinite(values) & values > 0);
if isempty(values)
    lim = 1;
    return;
end
maxValue = max(values);
values = sort(values(:));
idx = max(1, min(numel(values), round(percentile * numel(values))));
lim = min(maxValue, max(values(idx), floorFraction * maxValue));
if lim <= 0
    lim = 1;
end
end

function err = staggered_signed_error(coeff, kh)
offset = (0:numel(coeff)-1) + 0.5;
knumDelta = 2 * sin(kh * offset) * coeff(:);
err = knumDelta - kh;
end

function data = load_required_result(fileName)
if ~isfile(fileName)
    error('Missing required result file: %s', fileName);
end
data = load(fileName);
if ~isfield(data, 'result') || ~isfield(data.result, 'snapshots') || ~isfield(data.result.snapshots, 'vz')
    error('Result file has unexpected structure: %s', fileName);
end
end

function [standardCoeff, optimizedCoeff] = read_saved_coefficients(coeffFile, order)
if ~isfile(coeffFile)
    error('Coefficient file not found: %s', coeffFile);
end
fid = fopen(coeffFile, 'r', 'n', 'UTF-8');
cleanup = onCleanup(@() fclose(fid));
vectors = {};
while true
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end
    eq = strfind(line, '=');
    if isempty(eq)
        continue;
    end
    nums = sscanf(line(eq(1)+1:end), '%f').';
    if numel(nums) >= order
        vectors{end+1} = nums(1:order); %#ok<AGROW>
    end
end
if numel(vectors) < 2
    error('Could not read Taylor and optimized coefficient vectors from %s', coeffFile);
end
standardCoeff = vectors{1};
optimizedCoeff = vectors{2};
end

function write_summary_csv(outFile, summary, standardCoeff, optimizedCoeff)
fid = fopen(outFile, 'w', 'n', 'UTF-8');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '项目,数值\n');
fprintf(fid, '差分阶数,%d\n', summary.order);
fprintf(fid, 'khMax/pi,%.15g\n', summary.khMaxOverPi);
fprintf(fid, 'SA绝对误差阈值,%.15g\n', summary.targetError);
fprintf(fid, 'Taylor通带最大SA绝对误差,%.15g\n', summary.taylorMaxAbsError);
fprintf(fid, '优化系数通带最大SA绝对误差,%.15g\n', summary.optimizedMaxAbsError);
fprintf(fid, '改善倍数,%.15g\n', summary.improvement);
fprintf(fid, '误差定义,|k_num*Delta - k*Delta|\n');
fprintf(fid, '\n系数项,Taylor,优化系数\n');
for i = 1:numel(standardCoeff)
    fprintf(fid, 'c_%d,%.15g,%.15g\n', i, standardCoeff(i), optimizedCoeff(i));
end
end

function colors = standard_colors()
colors = struct();
colors.taylor = [0.86 0.25 0.10];
colors.optimized = [0.00 0.34 0.78];
end
