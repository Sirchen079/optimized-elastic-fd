% 沿从震源出发的固定方位角抽取波场空间测线，
% 对比 Taylor 系数与基于最大范数目标函数的优化系数：
%   (1) 主波峰振幅衰减；
%   (2) 主波峰相位（空间位置）延迟；
%   (3) 主波峰之后的拖尾振荡幅度。
% 仅读取已缓存的 .mat 结果，不重跑模拟、不修改主流程。

projectDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectDir, 'src'));
ewm_apply_chinese_style();

cfg = ewm_default_config(projectDir, 'standard');
cfg.output.figuresDir = fullfile(cfg.output.dir, 'figures');

model = ewm_load_marmousi(cfg.model);

stdData = load(fullfile(cfg.output.dir, 'exp3_staggered_pml_standard.mat'));
optData = load(fullfile(cfg.output.dir, 'exp3_staggered_pml_minimax.mat'));

% ------- 用户可调参数 -------
targetTime    = 1.250;   % 抽取测线的时刻 (s)
angleDeg      = 30;      % 测线相对水平方向的角度（>0 表示向下偏）
lineLengthM   = 6000;    % 测线从震源向外的长度 (m)
numSamples    = 1500;    % 沿测线的采样点数
trailWindowM  = 1500;    % 主波峰之后用于评估拖尾的窗长 (m)
trailGuardM   = 200;     % 主波峰附近保留的过渡带 (m)，避开主瓣污染
% ----------------------------

times = stdData.result.snapshots.time;
[~, snapIndex] = min(abs(times - targetTime));
actualTime = times(snapIndex);

taylorField    = stdData.result.snapshots.vz(:, :, snapIndex);
optimizedField = optData.result.snapshots.vz(:, :, snapIndex);

srcX = stdData.result.source.xM;
srcZ = stdData.result.source.depthM;

s = linspace(0, lineLengthM, numSamples);
xLine = srcX + s * cosd(angleDeg);
zLine = srcZ + s * sind(angleDeg);

inDomain = xLine >= model.x(1) & xLine <= model.x(end) ...
         & zLine >= model.z(1) & zLine <= model.z(end);
s = s(inDomain);
xLine = xLine(inDomain);
zLine = zLine(inDomain);

[Xgrid, Zgrid] = meshgrid(model.x, model.z);
taylorTrace    = interp2(Xgrid, Zgrid, taylorField,    xLine, zLine, 'linear', 0);
optimizedTrace = interp2(Xgrid, Zgrid, optimizedField, xLine, zLine, 'linear', 0);
diffTrace      = optimizedTrace - taylorTrace;

[absPeakT, idxPeakT] = max(abs(taylorTrace));
[absPeakO, idxPeakO] = max(abs(optimizedTrace));
peakSignedT = taylorTrace(idxPeakT);
peakSignedO = optimizedTrace(idxPeakO);
posPeakT = s(idxPeakT);
posPeakO = s(idxPeakO);
phaseDelayM = posPeakO - posPeakT;

% 拖尾区：主波峰之前（更接近震源）的一段，避开主瓣
trailLoT = posPeakT - trailGuardM - trailWindowM;
trailHiT = posPeakT - trailGuardM;
trailLoO = posPeakO - trailGuardM - trailWindowM;
trailHiO = posPeakO - trailGuardM;

trailMaskT = s >= trailLoT & s <= trailHiT;
trailMaskO = s >= trailLoO & s <= trailHiO;

if any(trailMaskT)
    trailRmsT  = sqrt(mean(taylorTrace(trailMaskT).^2));
    trailPeakT = max(abs(taylorTrace(trailMaskT)));
else
    trailRmsT = NaN; trailPeakT = NaN;
end
if any(trailMaskO)
    trailRmsO  = sqrt(mean(optimizedTrace(trailMaskO).^2));
    trailPeakO = max(abs(optimizedTrace(trailMaskO)));
else
    trailRmsO = NaN; trailPeakO = NaN;
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1400, 760]);
tl = tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
hT = plot(s / 1000, taylorTrace, 'Color', [0.10, 0.30, 0.85], ...
    'LineWidth', 1.4, 'DisplayName', 'Taylor 系数');
hold on;
hO = plot(s / 1000, optimizedTrace, '--', 'Color', [0.85, 0.20, 0.20], ...
    'LineWidth', 1.4, 'DisplayName', '最大范数目标函数优化系数');
yLim = get(gca, 'YLim');
plot([posPeakT, posPeakT] / 1000, yLim, ':', 'Color', [0.10, 0.30, 0.85], 'LineWidth', 0.9, 'HandleVisibility', 'off');
plot([posPeakO, posPeakO] / 1000, yLim, ':', 'Color', [0.85, 0.20, 0.20], 'LineWidth', 0.9, 'HandleVisibility', 'off');
xT = [trailLoT, trailHiT] / 1000;
patch([xT(1), xT(2), xT(2), xT(1)], [yLim(1), yLim(1), yLim(2), yLim(2)], ...
    [0.85, 0.85, 0.85], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
set(gca, 'YLim', yLim, 'Layer', 'top');
grid on; box on;
xlabel('沿测线距离 (km)');
ylabel('v_z (m/s)');
title(sprintf('沿震源方位 %g° 测线波场振幅，t = %.3f s', angleDeg, actualTime));
legend([hT, hO], 'Location', 'best');

nexttile;
plot(s / 1000, diffTrace, 'k-', 'LineWidth', 1.0);
grid on; box on;
xlabel('沿测线距离 (km)');
ylabel('差值 v_z (m/s)');
title('差值（优化 − Taylor）');

outFile = fullfile(cfg.output.figuresDir, ...
    sprintf('exp3_trace_angle%02ddeg_t%.3fs.png', round(angleDeg), actualTime));
ewm_save_figure(fig, outFile);
close(fig);

fprintf('已保存：%s\n', outFile);
fprintf('快照时间 = %.6f s（请求 %.3f s）\n', actualTime, targetTime);
fprintf('震源 (x, z) = (%.3f, %.3f) km，方位角 = %g°\n', srcX/1000, srcZ/1000, angleDeg);
fprintf('测线范围：s = 0–%.3f km（共 %d 个采样点）\n', s(end)/1000, numel(s));
fprintf('---- 主波峰 ----\n');
fprintf('  Taylor : |峰值| = %.4g m/s（带符号 %+.4g），位置 = %.4f km\n', ...
    absPeakT, peakSignedT, posPeakT/1000);
fprintf('  优化   : |峰值| = %.4g m/s（带符号 %+.4g），位置 = %.4f km\n', ...
    absPeakO, peakSignedO, posPeakO/1000);
fprintf('  振幅比 优化/Taylor = %.4f；相位延迟 = %+.4f km\n', ...
    absPeakO/absPeakT, phaseDelayM/1000);
fprintf('---- 拖尾区（主波峰前 %.2f–%.2f km）----\n', trailGuardM/1000, (trailGuardM+trailWindowM)/1000);
fprintf('  Taylor : RMS = %.4g m/s，峰值 = %.4g m/s\n', trailRmsT, trailPeakT);
fprintf('  优化   : RMS = %.4g m/s，峰值 = %.4g m/s\n', trailRmsO, trailPeakO);
fprintf('  比值 优化/Taylor: RMS = %.3f，峰值 = %.3f\n', ...
    trailRmsO/trailRmsT, trailPeakO/trailPeakT);
