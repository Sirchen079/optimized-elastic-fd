function ewm_plot_dispersion(standardCoeff, optimizedCoeff, metrics, cfg, optInfo)
%EWM_PLOT_DISPERSION 保存最大范数目标误差频散对比图。

if nargin < 5 || isempty(optInfo)
    optInfo = struct();
end

ewm_apply_chinese_style();
kh = metrics.kh(:);
x = kh / pi;
standardRatio = metrics.standard.ratio(:);
optimizedRatio = metrics.optimized.ratio(:);
standardError = standardRatio - 1;
optimizedError = optimizedRatio - 1;
standardAbsError = abs(standardError) .* kh;
optimizedAbsError = abs(optimizedError) .* kh;
target = cfg.coeff.targetError;

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1120, 860]);
tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(x, standardRatio, 'Color', [0.88 0.25 0.05], 'LineWidth', 1.8);
hold on;
plot(x, optimizedRatio, 'Color', [0.00 0.36 0.90], 'LineWidth', 1.8);
yline(1, 'k--', 'LineWidth', 1.2);
grid on;
xlim([0, cfg.coeff.khMax / pi]);
ylim([0.990, 1.002]);
xlabel('归一化波数 \theta / \pi');
ylabel('k_{num} / k');
title('数值波数曲线');
legend({'经典系数', '基于最大范数目标函数的优化系数', '精确值'}, ...
    'Location', 'south');

nexttile;
floorVal = max(target * 1.0e-2, 1.0e-7);
standardErrPlot = max(standardAbsError, floorVal);
optimizedErrPlot = max(optimizedAbsError, floorVal);
semilogy(x, standardErrPlot, 'Color', [0.88 0.25 0.05], 'LineWidth', 1.8);
hold on;
semilogy(x, optimizedErrPlot, 'Color', [0.00 0.36 0.90], 'LineWidth', 1.8);
yline(target, 'k--', 'LineWidth', 1.2);
grid on;
xlim([0, cfg.coeff.khMax / pi]);
ylim([floorVal, max([standardErrPlot; optimizedErrPlot]) * 2]);
xlabel('归一化波数 \theta / \pi');
ylabel('绝对误差 |k_{num}\Delta - k\Delta|');
title(sprintf('最大范数目标函数（绝对误差）：阈值 %.0e（log 轴）', target));
legend({'经典系数', '基于最大范数目标函数的优化系数', sprintf('%.0e 阈值', target)}, ...
    'Location', 'southeast');

mainTitle = sprintf('差分系数频散误差分析：kh ≤ %.2fπ，优化系数最大绝对误差 %.3g < %.0e（最大相对误差 %.3g）', ...
    cfg.coeff.khMax / pi, metrics.optimized.maxAbsError, target, metrics.optimized.maxRelError);
try
    sgtitle(fig, mainTitle, 'FontWeight', 'bold');
catch
    annotation(fig, 'textbox', [0.18, 0.955, 0.64, 0.035], ...
        'String', mainTitle, 'HorizontalAlignment', 'center', ...
        'EdgeColor', 'none', 'FontWeight', 'bold');
end

coeffFile = fullfile(cfg.output.dir, 'coefficients.txt');
fid = fopen(coeffFile, 'w');
fprintf(fid, '优化方法 = %s\n', optimizer_text(info_field(optInfo, 'method', 'unknown')));
fprintf(fid, '目标误差类型 = 绝对误差 max|k_num*Delta - k*Delta|（与 Zhang & Yao 2013 Eq. 23 一致）\n');
fprintf(fid, '目标误差 = %.15g\n', cfg.coeff.targetError);
fprintf(fid, '优化后最大绝对误差 = %.15g\n', metrics.optimized.maxAbsError);
fprintf(fid, '优化后最大相对误差 = %.15g\n', metrics.optimized.maxRelError);
fprintf(fid, '随机种子 = %.15g\n', info_field(optInfo, 'seed', NaN));
fprintf(fid, '优化采样点数 = %.15g\n', info_field(optInfo, 'samples', NaN));
fprintf(fid, '验证采样点数 = %.15g\n', info_field(optInfo, 'validationSamples', NaN));
fprintf(fid, 'Taylor系数 =');
fprintf(fid, ' %.15g', standardCoeff);
fprintf(fid, '\n优化系数 =');
fprintf(fid, ' %.15g', optimizedCoeff);
fprintf(fid, '\n优化系数逐项\n');
for m = 1:numel(optimizedCoeff)
    fprintf(fid, 'c_%d = %.15g\n', m, optimizedCoeff(m));
end
fprintf(fid, '\n最大kh = %.15g\n', cfg.coeff.khMax);
fprintf(fid, 'Taylor最大绝对误差 |k_num*Delta - k*Delta| = %.15g\n', metrics.standard.maxAbsError);
fprintf(fid, 'Taylor最大相对误差 |k_num/k - 1| = %.15g\n', metrics.standard.maxRelError);
fprintf(fid, '优化系数最大绝对误差 |k_num*Delta - k*Delta| = %.15g\n', metrics.optimized.maxAbsError);
fprintf(fid, '优化系数最大相对误差 |k_num/k - 1| = %.15g\n', metrics.optimized.maxRelError);
fprintf(fid, '绝对误差改善倍数 = %.15g\n', metrics.improvement);
fclose(fid);

ewm_save_figure(fig, fullfile(cfg.output.dir, 'figures', 'dispersion_curves.png'));
close(fig);
write_dispersion_csv(fullfile(cfg.output.dir, 'dispersion_curves.csv'), ...
    x, standardRatio, optimizedRatio, standardError, optimizedError);
end

function write_dispersion_csv(outFile, x, standardRatio, optimizedRatio, standardError, optimizedError)
fid = fopen(outFile, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '归一化波数theta除以pi,Taylor数值波数比,优化系数数值波数比,Taylor相对误差,优化系数相对误差\n');
for k = 1:numel(x)
    fprintf(fid, '%.15g,%.15g,%.15g,%.15g,%.15g\n', ...
        x(k), standardRatio(k), optimizedRatio(k), standardError(k), optimizedError(k));
end
end

function text = optimizer_text(value)
if isstring(value)
    value = char(value);
end
if strcmp(value, 'simulated_annealing_maximum_norm')
    text = '最大范数目标模拟退火';
elseif strcmp(value, 'unknown')
    text = '未知';
else
    text = value;
end
end

function value = info_field(info, fieldName, defaultValue)
if isstruct(info) && isfield(info, fieldName) && ~isempty(info.(fieldName))
    value = info.(fieldName);
else
    value = defaultValue;
end
end
