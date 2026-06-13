function metrics = ewm_dispersion_metrics(standardCoeff, optimizedCoeff, khMax, samples)
%EWM_DISPERSION_METRICS 对比导数算子的频散误差。
%
% 返回的 metrics 结构同时包含：
%   .maxAbsError —— 与 Zhang & Yao (2013) Eq. 23 一致的绝对误差最大范数
%                   max_{kh in (0, khMax]} | k_num*Delta - kh |，
%                   也是模拟退火优化所控制的目标量；
%   .maxRelError —— 数值波数比的最大偏差 max | k_num/k - 1 |，
%                   常用于交错网格的频散曲线绘制。

if nargin < 3 || isempty(khMax)
    khMax = 0.72 * pi;
end
if nargin < 4 || isempty(samples)
    samples = 1200;
end

kh = linspace(khMax / samples, khMax, samples).';

regularRatio = sin(kh) ./ kh;
standardRatio = staggered_ratio(standardCoeff, kh);
optimizedRatio = staggered_ratio(optimizedCoeff, kh);

metrics = struct();
metrics.kh = kh;
metrics.regular.ratio = regularRatio;
metrics.regular.maxRelError = max(abs(regularRatio - 1));
metrics.regular.maxAbsError = max(abs((regularRatio - 1) .* kh));
metrics.standard.ratio = standardRatio;
metrics.standard.maxRelError = max(abs(standardRatio - 1));
metrics.standard.maxAbsError = max(abs((standardRatio - 1) .* kh));
metrics.optimized.ratio = optimizedRatio;
metrics.optimized.maxRelError = max(abs(optimizedRatio - 1));
metrics.optimized.maxAbsError = max(abs((optimizedRatio - 1) .* kh));
metrics.improvement = metrics.standard.maxAbsError / metrics.optimized.maxAbsError;
end

function ratio = staggered_ratio(coeff, kh)
offset = (0:numel(coeff)-1) + 0.5;
ratio = 2 * sin(kh * offset) * coeff(:) ./ kh;
end
