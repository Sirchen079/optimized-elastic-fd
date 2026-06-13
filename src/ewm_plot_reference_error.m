function metrics = ewm_plot_reference_error(model, standardResult, minimaxResult, modelRef, referenceResult, outFile, metricsFile)
%EWM_PLOT_REFERENCE_ERROR 将两组波场与细网格参考解对比。

ewm_apply_chinese_style();

standard = standardResult.snapshots.vz;
minimax = minimaxResult.snapshots.vz;
reference = ewm_downsample_reference_snapshots(model, modelRef, referenceResult.snapshots.vz);

nSnap = min([size(standard, 3), size(minimax, 3), size(reference, 3)]);
times = standardResult.snapshots.time(1:nSnap);

stdErr = zeros(nSnap, 1);
minErr = zeros(nSnap, 1);
stdHighK = zeros(nSnap, 1);
minHighK = zeros(nSnap, 1);

figHeight = max(540, 360 * nSnap);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1500, figHeight]);
% 三张同色标的波场图：参考解、Taylor 残差、优化系数残差
layout = tiledlayout(fig, nSnap, 7, 'Padding', 'compact', 'TileSpacing', 'compact');

for k = 1:nSnap
    ref = reference(:, :, k);
    stdField = standard(:, :, k);
    minField = minimax(:, :, k);

    [stdErr(k), stdResidual] = normalized_reference_error(stdField, ref);
    [minErr(k), minResidual] = normalized_reference_error(minField, ref);
    stdHighK(k) = high_wavenumber_ratio(stdField);
    minHighK(k) = high_wavenumber_ratio(minField);

    % 用 98 分位数做色标截断，避免震源极值压制传播波前的对比度；
    % 三张波场图共用同一色标，让误差场以"参考波场偏淡版本"的形式呈现，
    % 视觉上正确反映出误差幅值小于波场幅值
    sharedLim = robust_clip([ref(:); stdField(:); minField(:)]);

    nexttile([1, 3]);
    imagesc(model.x / 1000, model.z / 1000, ref);
    axis image;
    set(gca, 'YDir', 'reverse');
    colormap(gca, ewm_wavefield_colormap());
    caxis([-sharedLim, sharedLim]);
    colorbar;
    title(sprintf('细网格参考，t = %.3f s', times(k)));
    xlabel('水平距离 (km)');
    ylabel('深度 (km)');

    nexttile([1, 2]);
    imagesc(model.x / 1000, model.z / 1000, stdResidual);
    axis image;
    set(gca, 'YDir', 'reverse');
    colormap(gca, ewm_wavefield_colormap());
    caxis([-sharedLim, sharedLim]);
    colorbar;
    title('Taylor 系数与参考解之差');
    xlabel('水平距离 (km)');
    ylabel('深度 (km)');

    nexttile([1, 2]);
    imagesc(model.x / 1000, model.z / 1000, minResidual);
    axis image;
    set(gca, 'YDir', 'reverse');
    colormap(gca, ewm_wavefield_colormap());
    caxis([-sharedLim, sharedLim]);
    colorbar;
    title('优化系数与参考解之差');
    xlabel('水平距离 (km)');
    ylabel('深度 (km)');
end

ewm_save_figure(fig, outFile);
close(fig);

metrics = struct();
metrics.time = times(:);
metrics.standardReferenceError = stdErr;
metrics.minimaxReferenceError = minErr;
metrics.standardHighKRatio = stdHighK;
metrics.minimaxHighKRatio = minHighK;
metrics.standardMeanReferenceError = mean(stdErr);
metrics.minimaxMeanReferenceError = mean(minErr);
metrics.standardMeanHighKRatio = mean(stdHighK);
metrics.minimaxMeanHighKRatio = mean(minHighK);
metrics.referenceSpacing = modelRef.dx;
metrics.comparisonSpacing = model.dx;
metrics.referenceNz = modelRef.nz;
metrics.referenceNx = modelRef.nx;
metrics.comparisonNz = model.nz;
metrics.comparisonNx = model.nx;

fid = fopen(metricsFile, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '参考解网格间距_m = %.15g\n', modelRef.dx);
fprintf(fid, '对比网格间距_m = %.15g\n', model.dx);
fprintf(fid, '说明：本表所列为细网格参考解与粗网格候选波场之间的相对 L2 偏离量，\n');
fprintf(fid, '      仅作为端到端波场可视化的辅助诊断；论文中报告的"误差"统一指\n');
fprintf(fid, '      模拟退火所控制的最大范数频散误差 max|k_num·Δ − k·Δ|（目标 1e-4）。\n\n');
fprintf(fid, '时间_s Taylor相对参考波场偏离 优化系数相对参考波场偏离 Taylor高波数能量占比 优化系数高波数能量占比\n');
for k = 1:nSnap
    fprintf(fid, '%.15g %.15g %.15g %.15g %.15g\n', ...
        times(k), stdErr(k), minErr(k), stdHighK(k), minHighK(k));
end
fprintf(fid, '\nTaylor平均相对参考波场偏离 = %.15g\n', metrics.standardMeanReferenceError);
fprintf(fid, '优化系数平均相对参考波场偏离 = %.15g\n', metrics.minimaxMeanReferenceError);
fprintf(fid, 'Taylor平均高波数能量占比 = %.15g\n', metrics.standardMeanHighKRatio);
fprintf(fid, '优化系数平均高波数能量占比 = %.15g\n', metrics.minimaxMeanHighKRatio);
end

function reference = ewm_downsample_reference_snapshots(model, modelRef, snapshots)
fz = round(model.dz / modelRef.dz);
fx = round(model.dx / modelRef.dx);
if fz < 1 || fx < 1
    error('参考解网格必须比对比网格更细。');
end

reference = snapshots(1:fz:end, 1:fx:end, :);
reference = reference(1:model.nz, 1:model.nx, :);
end

function [err, residual] = normalized_reference_error(candidate, reference)
den = sum(reference(:) .^ 2) + eps;
alpha = sum(candidate(:) .* reference(:)) / den;
scaledReference = alpha * reference;
residual = candidate - scaledReference;
err = norm(residual(:)) / (norm(scaledReference(:)) + eps);
end

function ratio = high_wavenumber_ratio(field)
field = field - mean(field(:));
spec = abs(fftshift(fft2(field))) .^ 2;
[nz, nx] = size(field);
kz = linspace(-1, 1, nz).';
kx = linspace(-1, 1, nx);
kr = sqrt(kz .^ 2 + kx .^ 2);
mask = kr > 0.55;
ratio = sum(spec(mask), 'all') / (sum(spec(:)) + eps);
end

function lim = robust_clip(values)
% 用 98 分位数代替 max(abs())，避免震源极值压制传播波前显示对比度。
values = abs(values(:));
values = values(isfinite(values) & values > 0);
if isempty(values)
    lim = 1;
    return;
end
maxValue = max(values);
sorted = sort(values);
idx = max(1, min(numel(sorted), round(0.98 * numel(sorted))));
percentileValue = sorted(idx);
% 不让百分位低于最大幅值的 18%，免得弱信号反过来被过度放大
lim = min(maxValue, max(percentileValue, 0.18 * maxValue));
if lim == 0
    lim = 1;
end
end

