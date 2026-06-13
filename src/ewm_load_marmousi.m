function model = ewm_load_marmousi(cfg)
%EWM_LOAD_MARMOUSI 加载 AGL 弹性 Marmousi 模型缓存。

spacing = cfg.spacing;
cacheFile = fullfile(cfg.root, sprintf('marmousi_cache_%dm.mat', spacing));

if exist(cacheFile, 'file')
    raw = load(cacheFile);
else
    baseFile = fullfile(cfg.root, 'marmousi_cache_20m.mat');
    if ~exist(baseFile, 'file')
        error('未找到 Marmousi 缓存文件，应存在：%s', baseFile);
    end
    raw = load(baseFile);
    if spacing < double(raw.dx)
        error('ewm:MissingFineMarmousiCache', ...
            '请求 %.3g m Marmousi 模型时必须存在对应缓存：%s。请先从原始 1.25 m SEG-Y 重采样生成该缓存。', ...
            spacing, cacheFile);
    end
    if abs(spacing / double(raw.dx) - round(spacing / double(raw.dx))) > 1.0e-9
        error('ewm:InvalidMarmousiSpacing', ...
            '请求的 Marmousi 间距 %.3g m 不是基础缓存 %.3g m 的整数倍。', spacing, double(raw.dx));
    end
    factor = max(1, round(spacing / double(raw.dx)));
    raw.vp = raw.vp(1:factor:end, 1:factor:end);
    raw.vs = raw.vs(1:factor:end, 1:factor:end);
    raw.rho = raw.rho(1:factor:end, 1:factor:end);
    raw.dx = double(raw.dx) * factor;
    raw.dz = double(raw.dz) * factor;
end

required = {'vp', 'vs', 'rho', 'dx', 'dz'};
for k = 1:numel(required)
    if ~isfield(raw, required{k})
        error('Marmousi 缓存缺少字段 "%s"。', required{k});
    end
end

vp = double(raw.vp);
vs = double(raw.vs);
rho = double(raw.rho);
dx = double(raw.dx);
dz = double(raw.dz);

if max(vp(:)) < 20
    vp = vp * 1000;
end
if max(vs(:)) > 0 && max(vs(:)) < 20
    vs = vs * 1000;
end
if median(rho(:), 'omitnan') < 20
    rho = rho * 1000;
end

vp(~isfinite(vp)) = median(vp(isfinite(vp)), 'omitnan');
vs(~isfinite(vs)) = 0;
rho(~isfinite(rho)) = median(rho(isfinite(rho)), 'omitnan');
vs = max(vs, 0);

if isfield(cfg, 'vsOffset') && cfg.vsOffset ~= 0
    vs = vs + cfg.vsOffset;
    vs = max(vs, 0);
end

if isfield(cfg, 'enforcePhysicalVs') && cfg.enforcePhysicalVs
    vsMax = 0.999 * vp / sqrt(2);
    vs = min(vs, vsMax);
end

if isfield(cfg, 'cropZ') && ~isempty(cfg.cropZ)
    vp = vp(cfg.cropZ, :);
    vs = vs(cfg.cropZ, :);
    rho = rho(cfg.cropZ, :);
end
if isfield(cfg, 'cropX') && ~isempty(cfg.cropX)
    vp = vp(:, cfg.cropX);
    vs = vs(:, cfg.cropX);
    rho = rho(:, cfg.cropX);
end

[nz, nx] = size(vp);

model = struct();
model.vp = vp;
model.vs = vs;
model.rho = rho;
model.dx = dx;
model.dz = dz;
model.nx = nx;
model.nz = nz;
model.x = (0:nx-1) * dx;
model.z = (0:nz-1) * dz;
if isfield(raw, 'sourceDescription') && ~isempty(raw.sourceDescription)
    model.source = raw.sourceDescription;
else
    model.source = 'AGL 弹性 Marmousi 缓存模型';
end
model.cacheSpacing = spacing;
if isfield(raw, 'sourceSpacing')
    model.sourceSpacing = raw.sourceSpacing;
else
    model.sourceSpacing = [];
end
if isfield(raw, 'resample')
    model.resample = raw.resample;
end
end
