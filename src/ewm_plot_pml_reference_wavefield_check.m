function ewm_plot_pml_reference_wavefield_check(model, noAbsorb, pml, referenceResult, cropInfo, outFile)
%EWM_PLOT_PML_REFERENCE_WAVEFIELD_CHECK 将 PML 波场与参考解对比。

ewm_apply_chinese_style();

ref = referenceResult.snapshots.vz(cropInfo.z, cropInfo.x, :);
noField = noAbsorb.snapshots.vz(:, :, end);
pmlField = pml.snapshots.vz(:, :, end);
refField = ref(:, :, end);
time = noAbsorb.snapshots.time(end);

noResidual = noField - refField;
pmlResidual = pmlField - refField;
noMinusPml = noField - pmlField;

waveLim = robust_limit([noField(:); pmlField(:); refField(:)], 0.995);
resLim = robust_limit([noResidual(:); pmlResidual(:); noMinusPml(:)], 0.995);

noErr = relative_l2(noResidual, refField);
pmlErr = relative_l2(pmlResidual, refField);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1520, 780]);
layout = tiledlayout(fig, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot_field(model, noField, waveLim, sprintf('（a）无吸收边界，t = %.3f s', time));

nexttile;
plot_field(model, pmlField, waveLim, sprintf('（b）PML 吸收边界，t = %.3f s', time));

nexttile;
plot_field(model, refField, waveLim, sprintf('（c）大域参考解，t = %.3f s', time));

nexttile;
plot_field(model, noResidual, resLim, sprintf('（d）无吸收 - 参考，E = %.3g', noErr));

nexttile;
plot_field(model, pmlResidual, resLim, sprintf('（e）PML - 参考，E = %.3g', pmlErr));

nexttile;
plot_field(model, noMinusPml, resLim, '（f）无吸收 - PML');

sgtitle(layout, 'PML 波场与大域参考解对照：用于判别底部波前是否为人工边界反射', ...
    'FontWeight', 'bold', 'FontSize', 15);

ewm_save_figure(fig, outFile);
close(fig);
end

function plot_field(model, field, limValue, panelTitle)
imagesc(model.x / 1000, model.z / 1000, field);
axis image;
set(gca, 'YDir', 'reverse');
colormap(gca, ewm_wavefield_colormap());
caxis([-limValue, limValue]);
colorbar;
title(panelTitle);
xlabel('水平距离 (km)');
ylabel('深度 (km)');
end

function limValue = robust_limit(values, fraction)
values = abs(values(:));
values = values(isfinite(values) & values > 0);
if isempty(values)
    limValue = 1;
    return;
end
values = sort(values);
idx = max(1, min(numel(values), round(fraction * numel(values))));
limValue = values(idx);
if limValue == 0
    limValue = max(values);
end
if limValue == 0
    limValue = 1;
end
end

function err = relative_l2(residual, reference)
err = norm(residual(:)) / (norm(reference(:)) + eps);
end
