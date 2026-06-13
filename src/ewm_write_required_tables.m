function tableFiles = ewm_write_required_tables(model, cfg, standardCoeff, optimizedCoeff, optInfo, dispersion, boundaryReferenceComparison, referenceComparison)
%EWM_WRITE_REQUIRED_TABLES 写入实验参数、参考解设置和 CFL 稳定性表格。

outDir = cfg.output.dir;

experimentRows = build_experiment_rows(model, cfg, standardCoeff, optimizedCoeff, optInfo, dispersion);
referenceRows = build_reference_rows(model, cfg, boundaryReferenceComparison, referenceComparison);
cflRows = build_cfl_rows(model, cfg, standardCoeff, optimizedCoeff, referenceComparison);

tableFiles = struct();
tableFiles.experimentParameters = write_table_artifacts(outDir, ...
    'experiment_parameters_table', {'参数', '取值'}, experimentRows, ...
    '完整实验参数表');
tableFiles.referenceSettings = write_table_artifacts(outDir, ...
    'reference_solution_settings_table', {'项目', '取值'}, referenceRows, ...
    '参考解设置表');
tableFiles.cflStability = write_table_artifacts(outDir, ...
    'cfl_stability_table', {'方案', '网格间距_m', '差分系数绝对值和', '稳定dt_s', '使用dt_s', '有效CFL', '裕度', '状态'}, cflRows, ...
    'CFL 稳定性表');
end

function rows = build_experiment_rows(model, cfg, standardCoeff, optimizedCoeff, optInfo, dispersion)
rows = {
    '运行模式', mode_text(cfg.mode);
    '模型来源', model.source;
    '模型缓存间距_m', model.cacheSpacing;
    '模型网格点数_z', model.nz;
    '模型网格点数_x', model.nx;
    '模型dz_m', model.dz;
    '模型dx_m', model.dx;
    '模型深度_km', max(model.z) / 1000;
    '模型宽度_km', max(model.x) / 1000;
    'Vp最小值_m_per_s', min(model.vp(:));
    'Vp最大值_m_per_s', max(model.vp(:));
    'Vs最小值_m_per_s', min(model.vs(:));
    'Vs最大值_m_per_s', max(model.vs(:));
    '密度最小值_kg_per_m3', min(model.rho(:));
    '密度最大值_kg_per_m3', max(model.rho(:));
    '差分半阶项数', cfg.coeff.order;
    '最大kh除以pi', cfg.coeff.khMax / pi;
    '误差定义', 'max|k_num·Δ − k·Δ| over kh ∈ [0, khMax]（与 SA 目标函数一致，全文唯一）';
    '频散误差目标_max_norm', cfg.coeff.targetError;
    'Taylor频散误差_max_norm', dispersion.standard.maxAbsError;
    '优化系数频散误差_max_norm', dispersion.optimized.maxAbsError;
    '频散误差改善倍数', dispersion.improvement;
    'CFL系数', cfg.sim.cfl;
    '时间步长_s', cfg.sim.dt;
    '时间步数', cfg.sim.nt;
    '总模拟时长_s', cfg.sim.totalTime;
    '最大快照时刻_s', cfg.sim.maxSnapshotTime;
    'PML层数', cfg.sim.nPml;
    'Ricker主频_Hz', cfg.sim.f0;
    'Ricker延迟周期数', cfg.sim.sourceDelayCycles;
    '震源振幅', cfg.sim.sourceAmplitude;
    '震源深度_m', cfg.sim.sourceDepthM;
    '震源水平位置比例', cfg.sim.sourceXFraction;
    '实验1快照时刻_s', join_numbers(cfg.snapshots.exp1Times);
    '实验2快照时刻_s', join_numbers(cfg.snapshots.exp2Times);
    '实验3快照时刻_s', join_numbers(cfg.snapshots.exp3Times);
    '优化方法', optimizer_text(info_field(optInfo, 'method', 'unknown'));
    '优化随机种子', info_field(optInfo, 'seed', NaN);
    '优化验证_max_norm误差', info_field(optInfo, 'validationMaxError', NaN);
    'Taylor系数', join_numbers(standardCoeff);
    '优化系数', join_numbers(optimizedCoeff)
    };
end

function rows = build_reference_rows(model, cfg, boundaryReferenceComparison, referenceComparison)
rows = {
    '对比模型网格点数_z', model.nz;
    '对比模型网格点数_x', model.nx;
    '对比模型间距_m', cfg.model.spacing;
    '实验2扩大计算域参考解启用', logical_text(cfg.boundaryReference.enabled);
    '实验2每侧扩展网格点数', cfg.boundaryReference.extraGrid;
    '实验2参考解PML层数', cfg.sim.nPml;
    '实验2参考解快照时刻_s', join_numbers(cfg.snapshots.exp2Times);
    '实验2参考解网格点数_z', value_or_field(boundaryReferenceComparison, 'referenceNz', NaN);
    '实验2参考解网格点数_x', value_or_field(boundaryReferenceComparison, 'referenceNx', NaN);
    '实验2波场偏离定义', '相对裁剪后的扩大计算域PML参考解的相对 L2 范数（辅助诊断，非论文"误差"）';
    '实验3细网格参考解启用', logical_text(cfg.reference.enabled);
    '实验3参考解目标间距_m', value_or_default(cfg.reference.spacing, NaN);
    '实验3对比模型间距_m', cfg.model.spacing;
    '实验3参考解快照时刻_s', join_numbers(cfg.snapshots.exp3Times);
    '实验3参考解实际间距_m', value_or_field(referenceComparison, 'referenceSpacing', NaN);
    '实验3参考解网格点数_z', value_or_field(referenceComparison, 'referenceNz', NaN);
    '实验3参考解网格点数_x', value_or_field(referenceComparison, 'referenceNx', NaN);
    '实验3参考解时间步长_s', value_or_field(referenceComparison, 'referenceDt', NaN);
    '实验3参考解时间步数', value_or_field(referenceComparison, 'referenceNt', NaN);
    '实验3参考解PML层数', value_or_field(referenceComparison, 'referencePmlLayers', NaN);
    '实验3参考解来源', value_or_field(referenceComparison, 'referenceSource', '未提供');
    '实验3原始数据间距_m', value_or_field(referenceComparison, 'referenceSourceSpacing', NaN);
    '实验3波场偏离定义', sprintf('相对%gm细网格参考解的相对 L2 范数（最小二乘幅值因子校正后），仅作辅助诊断，非论文"误差"', value_or_default(cfg.reference.spacing, NaN));
    '高波数能量占比定义', '径向波数大于0.55倍奈奎斯特波数的谱能量占总谱能量比例';
    '参考解空间间隔_m', value_or_default(cfg.reference.spacing, NaN);
    '参考解时间步长_s', value_or_field(referenceComparison, 'referenceDt', NaN);
    '参考解PML层数', value_or_field(referenceComparison, 'referencePmlLayers', NaN);
    '参考解网格规模', ref_grid_size_text(referenceComparison);
    '参考解是否重采样至粗网格', '是（空间上按最近邻下采样到主模型 nz × nx，时间上按最近时间步对齐）';
    '波场偏离统计区域是否排除PML', '否（仅裁剪到原物理域 nz × nx，未额外扣除 PML 厚度）';
    '比较分量', 'vz';
    '高波数阈值', '0.55 × Nyquist';
    '相对L2波场偏离是否进行最小二乘幅值校正', '是，先估计 alpha = <u, uref> / <uref, uref>，再计算 ||u - alpha uref||_2 / ||alpha uref||_2'
    };
end

function rows = build_cfl_rows(model, cfg, standardCoeff, optimizedCoeff, referenceComparison)
spacingLabel = sprintf('%gm主模型', model.dx);
rows = [
    cfl_row([spacingLabel, '_Taylor系数'], model.dx, model.dz, max(model.vp(:)), cfg.sim.cfl, cfg.sim.dt, standardCoeff);
    cfl_row([spacingLabel, '_优化系数'], model.dx, model.dz, max(model.vp(:)), cfg.sim.cfl, cfg.sim.dt, optimizedCoeff)
    ];

if isstruct(referenceComparison) && isfield(referenceComparison, 'referenceDt') && ...
        isfield(referenceComparison, 'referenceMaxVp') && isfield(referenceComparison, 'referenceSpacing')
    refSpacing = referenceComparison.referenceSpacing;
    refDt = referenceComparison.referenceDt;
    refLabel = sprintf('%gm参考解_Taylor系数', refSpacing);
    rows = [
        rows;
        cfl_row(refLabel, refSpacing, refSpacing, ...
            referenceComparison.referenceMaxVp, cfg.sim.cfl, refDt, standardCoeff)
        ];
end
end

function row = cfl_row(name, dx, dz, maxVp, cfl, usedDt, coeff)
dxMin = min(dx, dz);
gain = max(1, sum(abs(coeff)));
dtAtConfiguredCfl = cfl * dxMin / (sqrt(2) * maxVp * gain);
effectiveCfl = usedDt * sqrt(2) * maxVp * gain / dxMin;
margin = dtAtConfiguredCfl / (usedDt + eps);
if margin >= 1 - 1.0e-12
    status = '稳定';
else
    status = '不稳定';
end
row = {name, dxMin, gain, dtAtConfiguredCfl, usedDt, effectiveCfl, margin, status};
end

function files = write_table_artifacts(outDir, baseName, headers, rows, titleText)
figuresDir = fullfile(outDir, 'figures');
ewm_ensure_dir(figuresDir);
csvFile = fullfile(outDir, [baseName, '.csv']);
pngFile = fullfile(figuresDir, [baseName, '.png']);
write_csv(csvFile, headers, rows);
write_table_png(pngFile, headers, rows, titleText);
files = struct('csv', csvFile, 'png', pngFile);
end

function write_csv(outFile, headers, rows)
fid = fopen(outFile, 'w');
cleanup = onCleanup(@() fclose(fid));
write_csv_line(fid, headers);
for r = 1:size(rows, 1)
    row = cell(1, size(rows, 2));
    for c = 1:size(rows, 2)
        row{c} = value_to_text(rows{r, c});
    end
    write_csv_line(fid, row);
end
end

function write_csv_line(fid, row)
for c = 1:numel(row)
    if c > 1
        fprintf(fid, ',');
    end
    fprintf(fid, '%s', csv_escape(value_to_text(row{c})));
end
fprintf(fid, '\n');
end

function write_table_png(outFile, headers, rows, titleText)
ewm_apply_chinese_style();

rowCount = size(rows, 1);
figHeight = min(1500, max(360, 110 + 22 * (rowCount + 2)));
figWidth = 1500;
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, figWidth, figHeight]);
ax = axes('Parent', fig, 'Position', [0.025, 0.03, 0.95, 0.90]);
axis(ax, 'off');

text(ax, 0.5, 1.055, titleText, ...
    'Units', 'normalized', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', 15, ...
    'Interpreter', 'none');

lines = build_text_table_lines(headers, rows);
y = 0.99;
dy = min(0.042, 0.96 / max(1, numel(lines)));
for k = 1:numel(lines)
    text(ax, 0.01, y, lines{k}, ...
        'Units', 'normalized', 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 9.5, 'Interpreter', 'none');
    y = y - dy;
end

ewm_save_figure(fig, outFile);
close(fig);
end

function lines = build_text_table_lines(headers, rows)
textRows = cell(size(rows));
for r = 1:size(rows, 1)
    for c = 1:size(rows, 2)
        textRows{r, c} = value_to_text(rows{r, c});
    end
end

if size(rows, 2) == 2
    lines = cell(size(rows, 1) + 2, 1);
    lines{1} = sprintf('%-38s  %s', headers{1}, headers{2});
    lines{2} = [repmat('-', 1, 38), '  ', repmat('-', 1, 100)];
    for r = 1:size(rows, 1)
        lines{r + 2} = sprintf('%-38s  %s', textRows{r, 1}, textRows{r, 2});
    end
    return;
end

widths = [24, 12, 12, 14, 14, 12, 10, 8];
lines = cell(size(rows, 1) + 2, 1);
lines{1} = fixed_width_row(headers, widths);
lines{2} = fixed_width_row(repmat({'-'}, 1, numel(widths)), widths);
for r = 1:size(rows, 1)
    lines{r + 2} = fixed_width_row(textRows(r, :), widths);
end
end

function line = fixed_width_row(values, widths)
parts = cell(1, numel(widths));
for k = 1:numel(widths)
    text = value_to_text(values{k});
    if strcmp(text, '-')
        text = repmat('-', 1, widths(k));
    elseif numel(text) > widths(k)
        text = text(1:widths(k));
    end
    parts{k} = sprintf(['%-', num2str(widths(k)), 's'], text);
end
line = strjoin(parts, '  ');
end

function text = value_to_text(value)
if ischar(value)
    text = value;
elseif isstring(value)
    text = char(value);
elseif isnumeric(value)
    if isscalar(value)
        if isnan(value)
            text = '未提供';
        else
            text = sprintf('%.15g', value);
        end
    else
        text = join_numbers(value);
    end
elseif islogical(value)
    text = logical_text(value);
else
    text = '<不支持>';
end
end

function out = csv_escape(text)
needsQuote = contains(text, ',') || contains(text, '"') || contains(text, newline);
text = strrep(text, '"', '""');
if needsQuote
    out = ['"', text, '"'];
else
    out = text;
end
end

function out = join_numbers(values)
if isempty(values)
    out = '未提供';
    return;
end
values = values(:).';
parts = cell(1, numel(values));
for k = 1:numel(values)
    parts{k} = sprintf('%.15g', values(k));
end
out = strjoin(parts, ' ');
end

function out = logical_text(value)
if value
    out = '是';
else
    out = '否';
end
end

function value = info_field(info, fieldName, defaultValue)
if isstruct(info) && isfield(info, fieldName) && ~isempty(info.(fieldName))
    value = info.(fieldName);
else
    value = defaultValue;
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

function text = mode_text(value)
if isstring(value)
    value = char(value);
end
if strcmp(value, 'standard')
    text = '论文模式';
elseif strcmp(value, 'highfreq')
    text = '论文高频模式';
elseif strcmp(value, 'preview')
    text = '快速验证模式';
else
    text = value;
end
end

function value = value_or_field(s, fieldName, defaultValue)
if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end

function value = value_or_default(value, defaultValue)
if isempty(value)
    value = defaultValue;
end
end

function text = ref_grid_size_text(referenceComparison)
nz = value_or_field(referenceComparison, 'referenceNz', NaN);
nx = value_or_field(referenceComparison, 'referenceNx', NaN);
if isnan(nz) || isnan(nx)
    text = 'N/A';
else
    text = sprintf('%d x %d', nz, nx);
end
end
