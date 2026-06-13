function [coeff, info] = ewm_optimize_minimax_coeffs(order, khMax, samples, options)
%EWM_OPTIMIZE_MINIMAX_COEFFS 基于模拟退火的最大范数差分系数优化。
%
% 与 Zhang & Yao (2013) Eq. 23 一致：最小化绝对误差最大范数
%   E = max_{kh in [0, khMax]} | k*Delta - 2*sum_n c_n * sin((n+0.5)*kh) |
% 即对每个采样波数比较解析波数与差分算子数值波数的绝对偏差，
% 使用模拟退火替代最小二乘法求解该 minimax 问题。

if nargin < 1 || isempty(order)
    order = 5;
end
if nargin < 2 || isempty(khMax)
    khMax = 0.64 * pi;
end
if nargin < 3 || isempty(samples)
    samples = 1200;
end
if nargin < 4 || isempty(options)
    options = struct();
end

if order < 2
    error('ewm:InvalidOrder', '差分系数个数必须至少为 2。');
end

targetError = get_option(options, 'targetError', 1.0e-4);
if isfield(options, 'sa') && isstruct(options.sa)
    sa = options.sa;
else
    sa = struct();
end

standard = ewm_fd_coefficients('standard', order);
offset = (0:order-1) + 0.5;

kh = linspace(khMax / samples, khMax, samples).';
H = 2 * sin(kh * offset);
[base, reducedH] = reduced_symbol_matrix(H, offset);

seed = get_option(sa, 'seed', 20260425);
initialTemperature = get_option(sa, 'initialTemperature', 1.0);
minTemperature = get_option(sa, 'minTemperature', 1.0e-4);
coolingRate = get_option(sa, 'coolingRate', 0.992);
repeatPerTemperature = get_option(sa, 'repeatPerTemperature', 80);
restartCount = get_option(sa, 'restartCount', 20);
temperaturePower = get_option(sa, 'temperaturePower', 0.65);
acceptanceScale = get_option(sa, 'acceptanceScale', 8.0);
acceptanceFloor = get_option(sa, 'acceptanceFloor', 0.15);
minPerturbationFraction = get_option(sa, 'minPerturbationFraction', 0.002);
polishTrigger = get_option(sa, 'polishTrigger', 1.5);
validationSamples = get_option(sa, 'validationSamples', max(2000, samples));

defaultPerturbation = max(0.35 * abs(standard(2:end)), 0.02 ./ (1:order-1) .^ 2);
perturbationScale = row_option(sa, 'perturbationScale', defaultPerturbation, order - 1);
restartScale = row_option(sa, 'restartScale', 0.7 * perturbationScale, order - 1);

standardPack = standard(2:end);
bestPack = standardPack;
bestObjective = maximum_error(bestPack, base, reducedH, offset, kh);
startObjective = bestObjective;

oldRng = rng;
cleanup = onCleanup(@() rng(oldRng));
rng(seed, 'twister');

evalCount = 0;
temperatureSteps = 0;
polishEvals = 0;

traceCapacity = 256;
traceRestart = zeros(traceCapacity, 1);
traceStep = zeros(traceCapacity, 1);
traceTemperature = zeros(traceCapacity, 1);
traceCurrentObjective = zeros(traceCapacity, 1);
traceBestObjective = zeros(traceCapacity, 1);
traceCount = 0;
traceRestart(1) = 1;
traceStep(1) = 0;
traceTemperature(1) = initialTemperature;
traceCurrentObjective(1) = bestObjective;
traceBestObjective(1) = bestObjective;
traceCount = 1;

for restart = 1:restartCount
    if restart == 1
        currentPack = standardPack;
    else
        currentPack = random_valid_pack(bestPack, restartScale, offset);
    end
    currentObjective = maximum_error(currentPack, base, reducedH, offset, kh);

    temperature = initialTemperature;
    while temperature > minTemperature
        temperatureRatio = max(temperature / initialTemperature, minTemperature / initialTemperature);
        step = perturbationScale .* temperatureRatio .^ temperaturePower + ...
            perturbationScale .* minPerturbationFraction;
        acceptTemperature = targetError * acceptanceScale * (acceptanceFloor + temperatureRatio);

        for n = 1:repeatPerTemperature
            candidatePack = currentPack + randn(1, order - 1) .* step;
            if ~valid_coefficients(unpack_coeff(candidatePack, offset))
                continue;
            end

            candidateObjective = maximum_error(candidatePack, base, reducedH, offset, kh);
            evalCount = evalCount + 1;
            if candidateObjective < currentObjective || ...
                    rand() < exp((currentObjective - candidateObjective) / max(acceptTemperature, realmin))
                currentPack = candidatePack;
                currentObjective = candidateObjective;

                if currentObjective < bestObjective
                    bestPack = currentPack;
                    bestObjective = currentObjective;
                end
            end
        end

        temperature = temperature * coolingRate;
        temperatureSteps = temperatureSteps + 1;

        traceCount = traceCount + 1;
        if traceCount > traceCapacity
            traceCapacity = traceCapacity * 2;
            traceRestart(traceCapacity) = 0;
            traceStep(traceCapacity) = 0;
            traceTemperature(traceCapacity) = 0;
            traceCurrentObjective(traceCapacity) = 0;
            traceBestObjective(traceCapacity) = 0;
        end
        traceRestart(traceCount) = restart;
        traceStep(traceCount) = temperatureSteps;
        traceTemperature(traceCount) = temperature;
        traceCurrentObjective(traceCount) = currentObjective;
        traceBestObjective(traceCount) = bestObjective;
    end

    if bestObjective <= polishTrigger * targetError
        [bestPack, bestObjective, usedEvals] = polish_pack(bestPack, bestObjective, ...
            base, reducedH, offset, kh, sa);
        polishEvals = polishEvals + usedEvals;

        traceCount = traceCount + 1;
        if traceCount > traceCapacity
            traceCapacity = traceCapacity * 2;
            traceRestart(traceCapacity) = 0;
            traceStep(traceCapacity) = 0;
            traceTemperature(traceCapacity) = 0;
            traceCurrentObjective(traceCapacity) = 0;
            traceBestObjective(traceCapacity) = 0;
        end
        traceRestart(traceCount) = restart;
        traceStep(traceCount) = temperatureSteps;
        traceTemperature(traceCount) = -1;
        traceCurrentObjective(traceCount) = bestObjective;
        traceBestObjective(traceCount) = bestObjective;
    end

    if bestObjective <= targetError
        break;
    end
end

coeff = unpack_coeff(bestPack, offset);

validationKh = linspace(khMax / validationSamples, khMax, validationSamples).';
validationH = 2 * sin(validationKh * offset);
[validationBase, validationReducedH] = reduced_symbol_matrix(validationH, offset);
validationObjective = maximum_error(bestPack, validationBase, validationReducedH, offset, validationKh);
if validationObjective > targetError
    [bestPack, validationObjective, usedEvals] = polish_pack(bestPack, validationObjective, ...
        validationBase, validationReducedH, offset, validationKh, sa);
    polishEvals = polishEvals + usedEvals;
    coeff = unpack_coeff(bestPack, offset);
    bestObjective = maximum_error(bestPack, base, reducedH, offset, kh);
end

if validationObjective > targetError
    error('ewm:AnnealingTargetNotMet', ...
        ['模拟退火法未达到目标误差 %.3g；当前最大误差 %.6g。', ...
         '请增加 restartCount/repeatPerTemperature，或减小 khMax。'], ...
        targetError, validationObjective);
end

info = struct();
info.method = 'simulated_annealing_maximum_norm';
info.order = order;
info.khMax = khMax;
info.samples = samples;
info.validationSamples = validationSamples;
info.targetError = targetError;
info.seed = seed;
info.startMaxError = startObjective;
info.objective = bestObjective;
info.validationMaxError = validationObjective;
info.evalCount = evalCount;
info.polishEvalCount = polishEvals;
info.temperatureSteps = temperatureSteps;
info.restartCount = restartCount;
info.initialTemperature = initialTemperature;
info.minTemperature = minTemperature;
info.coolingRate = coolingRate;
info.repeatPerTemperature = repeatPerTemperature;

trace = struct();
trace.restart = traceRestart(1:traceCount);
trace.step = traceStep(1:traceCount);
trace.temperature = traceTemperature(1:traceCount);
trace.currentObjective = traceCurrentObjective(1:traceCount);
trace.bestObjective = traceBestObjective(1:traceCount);
info.trace = trace;
end

function [base, reducedH] = reduced_symbol_matrix(H, offset)
% H 的列对应 2*sin(kh*(n+0.5))（不再除以 kh），因此
% base + reducedH*pack(:) 直接给出离散一阶导数算子在 kh 处的数值波数 k_num*Delta。
% 低波数一致性约束通过把 c1 表达为 c2..c_{N/2} 的线性组合（保证 2*sum(offset.*c)=1）实现，
% 这一关系在 unpack_coeff 中显式写出，此处对 H 做相同的代数消元。
c1Base = 1 / (2 * offset(1));
c1Weights = -offset(2:end) / offset(1);
base = H(:, 1) * c1Base;
reducedH = H(:, 2:end) + H(:, 1) * c1Weights;
end

function val = maximum_error(pack, base, reducedH, offset, target)
% 与 Zhang & Yao (2013) Eq. 23 一致：返回绝对误差的最大范数
%   max_{kh} | k_num*Delta - kh |.
coeff = unpack_coeff(pack, offset);
if ~valid_coefficients(coeff)
    val = inf;
    return;
end

symbol = base + reducedH * pack(:);
val = max(abs(symbol - target(:)));
end

function coeff = unpack_coeff(pack, offset)
coeff = zeros(1, numel(offset));
coeff(2:end) = pack(:).';
coeff(1) = (1 - 2 * sum(offset(2:end) .* coeff(2:end))) / (2 * offset(1));
end

function tf = valid_coefficients(coeff)
signPattern = (-1) .^ (0:numel(coeff)-1);
if any(coeff .* signPattern <= 0)
    tf = false;
    return;
end

absCoeff = abs(coeff);
tf = all(absCoeff(1:end-1) > absCoeff(2:end));
end

function pack = random_valid_pack(centerPack, scale, offset)
for trial = 1:2000
    pack = centerPack + randn(size(centerPack)) .* scale;
    if valid_coefficients(unpack_coeff(pack, offset))
        return;
    end
end

pack = centerPack;
end

function [bestPack, bestObjective, evalCount] = polish_pack(bestPack, bestObjective, ...
    base, reducedH, offset, target, sa)
evalCount = 0;
nvar = numel(bestPack);
scaleList = row_option(sa, 'polishScales', ...
    [1.0e-3, 5.0e-4, 2.0e-4, 1.0e-4, 5.0e-5, 2.0e-5, ...
     1.0e-5, 5.0e-6, 2.0e-6, 1.0e-6, 5.0e-7, 2.0e-7, 1.0e-7], 13);
stepRatios = row_option(sa, 'polishStepRatios', 1 ./ (1:nvar) .^ 1.8, nvar);
trialsPerScale = get_option(sa, 'polishTrialsPerScale', 2500);
maxPasses = get_option(sa, 'polishMaxPasses', 4);

for s = scaleList
    improved = true;
    pass = 0;
    while improved && pass < maxPasses
        improved = false;
        pass = pass + 1;
        step = s .* stepRatios;

        for n = 1:trialsPerScale
            candidatePack = bestPack + randn(1, nvar) .* step;
            if ~valid_coefficients(unpack_coeff(candidatePack, offset))
                continue;
            end

            candidateObjective = maximum_error(candidatePack, base, reducedH, offset, target);
            evalCount = evalCount + 1;
            if candidateObjective < bestObjective
                bestPack = candidatePack;
                bestObjective = candidateObjective;
                improved = true;
            end
        end
    end
end
end

function value = get_option(options, fieldName, defaultValue)
if isstruct(options) && isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end

function value = row_option(options, fieldName, defaultValue, expectedLength)
value = get_option(options, fieldName, defaultValue);
value = value(:).';
if numel(value) == 1 && expectedLength > 1
    value = repmat(value, 1, expectedLength);
elseif numel(value) ~= expectedLength
    error('ewm:InvalidAnnealingOption', ...
        '模拟退火参数 %s 的长度应为 %d。', fieldName, expectedLength);
end
end
