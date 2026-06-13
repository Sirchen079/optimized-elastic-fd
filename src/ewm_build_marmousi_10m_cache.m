function cacheFile = ewm_build_marmousi_10m_cache(modelRoot, force)
%EWM_BUILD_MARMOUSI_10M_CACHE 从原始 1.25 m SEG-Y 构建 10 m Marmousi 缓存。

if nargin < 2 || isempty(force)
    force = false;
end

targetSpacing = 10;
sourceSpacing = 1.25;
decimation = round(targetSpacing / sourceSpacing);
if abs(decimation * sourceSpacing - targetSpacing) > 1.0e-9
    error('ewm:InvalidResampling', '10 m 必须是原始网格间距的整数倍。');
end

cacheFile = fullfile(modelRoot, 'marmousi_cache_10m.mat');
sourceFiles = struct();
sourceFiles.vp = fullfile(modelRoot, 'MODEL_P-WAVE_VELOCITY_1.25m.segy');
sourceFiles.vs = fullfile(modelRoot, 'MODEL_S-WAVE_VELOCITY_1.25m.segy');
sourceFiles.rho = fullfile(modelRoot, 'MODEL_DENSITY_1.25m.segy');

assert_source_exists(sourceFiles.vp);
assert_source_exists(sourceFiles.vs);
assert_source_exists(sourceFiles.rho);

if ~force && cache_is_current(cacheFile, sourceFiles)
    return;
end

fprintf('正在从原始 1.25 m SEG-Y 重采样生成 10 m Marmousi 缓存...\n');
fprintf('  读取 Vp：%s\n', sourceFiles.vp);
[vp, vpInfo] = read_decimated_segy(sourceFiles.vp, decimation);
fprintf('  读取 Vs：%s\n', sourceFiles.vs);
[vs, vsInfo] = read_decimated_segy(sourceFiles.vs, decimation);
fprintf('  读取 rho：%s\n', sourceFiles.rho);
[rho, rhoInfo] = read_decimated_segy(sourceFiles.rho, decimation);

if ~isequal(size(vp), size(vs), size(rho))
    error('ewm:MarmousiGridMismatch', '重采样后的 Vp、Vs、rho 网格尺寸不一致。');
end
if vpInfo.nSamples ~= vsInfo.nSamples || vpInfo.nSamples ~= rhoInfo.nSamples || ...
        vpInfo.nTraces ~= vsInfo.nTraces || vpInfo.nTraces ~= rhoInfo.nTraces
    error('ewm:MarmousiSegyMismatch', '三个原始 SEG-Y 文件的网格布局不一致。');
end

dx = targetSpacing;
dz = targetSpacing;
resample = struct();
resample.method = 'raw_1p25m_segy_decimation';
resample.sourceSpacing = sourceSpacing;
resample.targetSpacing = targetSpacing;
resample.decimation = decimation;
resample.sourceSamplesPerTrace = vpInfo.nSamples;
resample.sourceTraceCount = vpInfo.nTraces;
resample.outputNz = size(vp, 1);
resample.outputNx = size(vp, 2);
resample.formatCode = vpInfo.formatCode;
resample.sampleInterval = vpInfo.sampleInterval;
resample.sourceFiles = sourceFiles;
sourceDescription = '原始 1.25 m SEG-Y 重采样得到的 10 m AGL 弹性 Marmousi 模型';

fprintf('  10 m 缓存尺寸：nz=%d, nx=%d\n', size(vp, 1), size(vp, 2));
save(cacheFile, 'vp', 'vs', 'rho', 'dx', 'dz', 'sourceSpacing', ...
    'sourceFiles', 'sourceDescription', 'resample', '-v7.3');
fprintf('  已写入：%s\n', cacheFile);
end

function assert_source_exists(fileName)
if ~exist(fileName, 'file')
    error('ewm:MissingMarmousiSegy', '缺少原始 SEG-Y 文件：%s', fileName);
end
end

function ok = cache_is_current(cacheFile, sourceFiles)
ok = false;
if ~exist(cacheFile, 'file')
    return;
end

cacheInfo = dir(cacheFile);
sources = {sourceFiles.vp, sourceFiles.vs, sourceFiles.rho};
for k = 1:numel(sources)
    info = dir(sources{k});
    if isempty(info) || info.datenum > cacheInfo.datenum
        return;
    end
end

try
    raw = load(cacheFile, 'dx', 'dz', 'vp', 'vs', 'rho', 'resample');
    ok = isequal(size(raw.vp), size(raw.vs), size(raw.rho)) && ...
        abs(double(raw.dx) - 10) < eps && abs(double(raw.dz) - 10) < eps;
    if ok && isfield(raw, 'resample') && isfield(raw.resample, 'method')
        ok = strcmp(raw.resample.method, 'raw_1p25m_segy_decimation');
    end
catch
    ok = false;
end
end

function [field, info] = read_decimated_segy(fileName, decimation)
info = inspect_segy(fileName);
if info.formatCode ~= 1 && info.formatCode ~= 5
    error('ewm:UnsupportedSegyFormat', ...
        '只支持 SEG-Y 格式码 1(IBM float) 或 5(IEEE float)，当前文件格式码为 %d：%s', ...
        info.formatCode, fileName);
end

sampleIdx = 1:decimation:info.nSamples;
traceIdx = 1:decimation:info.nTraces;
field = zeros(numel(sampleIdx), numel(traceIdx));

fid = fopen(fileName, 'r', 'ieee-be');
if fid < 0
    error('ewm:FileOpenFailed', '无法打开 SEG-Y 文件：%s', fileName);
end
cleanup = onCleanup(@() fclose(fid));

for j = 1:numel(traceIdx)
    traceNumber = traceIdx(j);
    traceOffset = 3600 + (traceNumber - 1) * info.traceBytes;
    status = fseek(fid, traceOffset + 240, 'bof');
    if status ~= 0
        error('ewm:SegySeekFailed', '无法定位 SEG-Y 道数据：%s，道号 %d', fileName, traceNumber);
    end

    if info.formatCode == 1
        bytes = fread(fid, [4, info.nSamples], 'uint8=>uint8');
        if size(bytes, 2) ~= info.nSamples
            error('ewm:SegyReadFailed', 'SEG-Y 道数据读取不完整：%s，道号 %d', fileName, traceNumber);
        end
        trace = ibm_float32_to_double(bytes);
    else
        trace = fread(fid, info.nSamples, 'single=>double');
        if numel(trace) ~= info.nSamples
            error('ewm:SegyReadFailed', 'SEG-Y 道数据读取不完整：%s，道号 %d', fileName, traceNumber);
        end
        trace = trace(:).';
    end

    field(:, j) = trace(sampleIdx).';
end
end

function info = inspect_segy(fileName)
fileInfo = dir(fileName);
if isempty(fileInfo)
    error('ewm:MissingMarmousiSegy', '缺少原始 SEG-Y 文件：%s', fileName);
end

fid = fopen(fileName, 'r', 'ieee-be');
if fid < 0
    error('ewm:FileOpenFailed', '无法打开 SEG-Y 文件：%s', fileName);
end
cleanup = onCleanup(@() fclose(fid));

fseek(fid, 3200 + 16, 'bof');
sampleInterval = fread(fid, 1, 'uint16=>double');
fseek(fid, 3200 + 20, 'bof');
nSamples = fread(fid, 1, 'uint16=>double');
fseek(fid, 3200 + 24, 'bof');
formatCode = fread(fid, 1, 'uint16=>double');

bytesPerSample = 4;
traceBytes = 240 + bytesPerSample * nSamples;
dataBytes = fileInfo.bytes - 3600;
nTraces = dataBytes / traceBytes;
if nSamples <= 0 || abs(nTraces - round(nTraces)) > 0
    error('ewm:InvalidSegyLayout', '无法识别 SEG-Y 网格布局：%s', fileName);
end

info = struct();
info.sampleInterval = sampleInterval;
info.nSamples = nSamples;
info.formatCode = formatCode;
info.bytesPerSample = bytesPerSample;
info.traceBytes = traceBytes;
info.nTraces = round(nTraces);
end

function values = ibm_float32_to_double(bytes)
first = bytes(1, :);
negative = bitand(first, uint8(128)) ~= 0;
exponent = double(bitand(first, uint8(127))) - 64;
mantissa = double(bytes(2, :)) * 65536 + double(bytes(3, :)) * 256 + double(bytes(4, :));

values = (mantissa / 16777216) .* 16 .^ exponent;
values(negative) = -values(negative);
values(mantissa == 0) = 0;
end
