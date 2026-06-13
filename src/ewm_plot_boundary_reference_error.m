function metrics = ewm_plot_boundary_reference_error(model, noAbsorb, pml, modelRef, referenceResult, cropInfo, outFile, metricsFile)
%EWM_PLOT_BOUNDARY_REFERENCE_ERROR 将边界方案与扩大计算域参考解对比。

ewm_apply_chinese_style();

ref = referenceResult.snapshots.vz(cropInfo.z, cropInfo.x, :);
noAbsorbField = noAbsorb.snapshots.vz;
pmlField = pml.snapshots.vz;
nSnap = min([size(ref, 3), size(noAbsorbField, 3), size(pmlField, 3)]);
times = noAbsorb.snapshots.time(1:nSnap);

noAbsErr = zeros(nSnap, 1);
pmlErr = zeros(nSnap, 1);
noAbsBoundaryErr = zeros(nSnap, 1);
pmlBoundaryErr = zeros(nSnap, 1);

figHeight = max(720, 520 * nSnap);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1520, figHeight]);
layout = tiledlayout(fig, 2 * nSnap, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for k = 1:nSnap
    refK = ref(:, :, k);
    noK = noAbsorbField(:, :, k);
    pmlK = pmlField(:, :, k);
    noRes = noK - refK;
    pmlRes = pmlK - refK;

    noAbsErr(k) = relative_l2(noRes, refK);
    pmlErr(k) = relative_l2(pmlRes, refK);
    noAbsBoundaryErr(k) = relative_l2_boundary(noRes, refK);
    pmlBoundaryErr(k) = relative_l2_boundary(pmlRes, refK);

    waveLim = max(abs([refK(:); noK(:); pmlK(:)]));
    errLim = max(abs([noRes(:); pmlRes(:)]));
    if waveLim == 0
        waveLim = 1;
    end
    if errLim == 0
        errLim = 1;
    end

    nexttile;
    plot_field(model, refK, waveLim, sprintf('（a）大域参考解，t = %.3f s', times(k)));

    nexttile;
    plot_field(model, noRes, errLim, '（b）无吸收边界与参考解之差');

    nexttile;
    plot_field(model, pmlRes, errLim, '（c）PML 边界与参考解之差');

    nexttile([1, 3]);
    metricValues = [noAbsErr(k), pmlErr(k); noAbsBoundaryErr(k), pmlBoundaryErr(k)];
    b = bar(metricValues, 0.62);
    b(1).FaceColor = [0.05, 0.32, 0.62];
    b(2).FaceColor = [0.82, 0.28, 0.12];
    set(gca, 'XTickLabel', {'全域', '边界带'});
    style_axes();
    ylabel('相对 L2 波场偏离');
    reductionFull = percent_reduction(noAbsErr(k), pmlErr(k));
    reductionBoundary = percent_reduction(noAbsBoundaryErr(k), pmlBoundaryErr(k));
    title(sprintf('（d）相对参考解的波场偏离：PML 全域降幅 %.1f%%，边界带降幅 %.1f%%', ...
        reductionFull, reductionBoundary));
    legend({'无吸收边界', 'PML 吸收边界'}, ...
        'Location', 'northoutside', 'Orientation', 'horizontal');
    add_grouped_bar_values(metricValues);
    ylim([0, max(metricValues(:)) * 1.18 + eps]);
end

ewm_save_figure(fig, outFile);
close(fig);

metrics = struct();
metrics.time = times(:);
metrics.noAbsorbReferenceError = noAbsErr;
metrics.pmlReferenceError = pmlErr;
metrics.noAbsorbBoundaryReferenceError = noAbsBoundaryErr;
metrics.pmlBoundaryReferenceError = pmlBoundaryErr;
metrics.noAbsorbMeanReferenceError = mean(noAbsErr);
metrics.pmlMeanReferenceError = mean(pmlErr);
metrics.noAbsorbMeanBoundaryReferenceError = mean(noAbsBoundaryErr);
metrics.pmlMeanBoundaryReferenceError = mean(pmlBoundaryErr);
metrics.extraGrid = cropInfo.extraGrid;
metrics.referenceNx = modelRef.nx;
metrics.referenceNz = modelRef.nz;

fid = fopen(metricsFile, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '参考解类型 = 扩大计算域PML参考解，并裁剪到原始模型范围\n');
fprintf(fid, '每侧扩展网格点数 = %d\n', cropInfo.extraGrid);
fprintf(fid, '参考解网格点数_z = %d\n', modelRef.nz);
fprintf(fid, '参考解网格点数_x = %d\n', modelRef.nx);
fprintf(fid, '说明：表中数值为候选波场相对大域参考解的相对 L2 偏离量，仅用于评估\n');
fprintf(fid, '      PML 与无吸收边界的反射差异；论文中报告的"误差"统一指模拟退火\n');
fprintf(fid, '      所控制的最大范数频散误差 max|k_num·Δ − k·Δ|（目标 1e-4）。\n\n');
fprintf(fid, '时间_s 无吸收全域偏离 PML全域偏离 无吸收边界带偏离 PML边界带偏离\n');
for k = 1:nSnap
    fprintf(fid, '%.15g %.15g %.15g %.15g %.15g\n', ...
        times(k), noAbsErr(k), pmlErr(k), noAbsBoundaryErr(k), pmlBoundaryErr(k));
end
fprintf(fid, '\n无吸收平均全域波场偏离 = %.15g\n', metrics.noAbsorbMeanReferenceError);
fprintf(fid, 'PML平均全域波场偏离 = %.15g\n', metrics.pmlMeanReferenceError);
fprintf(fid, '无吸收平均边界带波场偏离 = %.15g\n', metrics.noAbsorbMeanBoundaryReferenceError);
fprintf(fid, 'PML平均边界带波场偏离 = %.15g\n', metrics.pmlMeanBoundaryReferenceError);
end

function plot_field(model, field, limValue, panelTitle)
imagesc(model.x / 1000, model.z / 1000, field);
axis tight;
set(gca, 'YDir', 'reverse');
colormap(gca, ewm_wavefield_colormap());
caxis([-limValue, limValue]);
colorbar;
style_axes();
title(panelTitle);
xlabel('水平距离 (km)');
ylabel('深度 (km)');
end

function err = relative_l2(residual, reference)
err = norm(residual(:)) / (norm(reference(:)) + eps);
end

function err = relative_l2_boundary(residual, reference)
[nz, nx] = size(residual);
width = max(3, round(0.18 * min(nz, nx)));
mask = false(nz, nx);
mask(1:width, :) = true;
mask(end-width+1:end, :) = true;
mask(:, 1:width) = true;
mask(:, end-width+1:end) = true;
err = norm(residual(mask)) / (norm(reference(mask)) + eps);
end

function reduction = percent_reduction(baseline, candidate)
reduction = 100 * (baseline - candidate) / (baseline + eps);
end

function add_grouped_bar_values(values)
[nGroup, nSeries] = size(values);
groupWidth = min(0.8, nSeries / (nSeries + 1.5));
for i = 1:nGroup
    for j = 1:nSeries
        x = i - groupWidth / 2 + (2 * j - 1) * groupWidth / (2 * nSeries);
        text(x, values(i, j), sprintf('%.3g', values(i, j)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontWeight', 'bold', ...
            'FontSize', 10);
    end
end
end

function style_axes()
grid on;
box on;
set(gca, 'LineWidth', 0.9, 'FontSize', 11, 'GridAlpha', 0.18);
end
