% 一维高斯脉冲传播对比图：
%   将 Taylor 系数与优化系数（最大范数目标函数）得到的有限差分解
%   分别与解析解（即初始波形的精确平移）对比，
%   以单点扫描方式重现“振幅 vs x/Δx”的经典色散评估图。
%
% 该脚本独立运行：只调用 src/ 中的现有系数工具，不修改主流程，
% 也不依赖二维 Marmousi 缓存。优化系数从缓存 .mat 读取 result.coeff。

projectDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectDir, 'src'));
ewm_apply_chinese_style();

% ------- 参数 -------
order      = 5;         % 半阶数，对应 10 阶空间精度
N          = 320;       % 网格点数
dx         = 1.0;       % 网格间距（无量纲）
c          = 1.0;       % 波速（无量纲）
cflFactor  = 0.1;       % 时间步安全系数（设小以隔离空间色散，避免时间误差掩盖差异）
x0         = 10;        % 初始脉冲中心位置（以 dx 为单位）
sigma      = 1.5;       % 高斯脉冲半宽（以 dx 为单位）；窄脉冲含高波数成分
amp        = 0.5;       % 脉冲幅值
finalPos   = 220;       % 最终脉冲中心位置（以 dx 为单位）
viewWin    = [200, 230];% 绘图窗口（x/Δx）
% --------------------

taylorCoeff = ewm_fd_coefficients('standard', order);

minimaxFile = fullfile(projectDir, 'results_standard', 'exp3_staggered_pml_minimax.mat');
optData = load(minimaxFile);
optCoeff = optData.result.coeff(:).';
assert(numel(optCoeff) == order, '优化系数阶数与脚本设定不一致。');

dtTaylor = cflFactor * dx / (c * sum(abs(taylorCoeff)));
dtOpt    = cflFactor * dx / (c * sum(abs(optCoeff)));
dt = min(dtTaylor, dtOpt);  % 用同一 dt，使两条结果对齐

totalDist = (finalPos - x0) * dx;
nt = round(totalDist / (c * dt));
elapsed = nt * dt;

x = (0:N-1) * dx;            % p 的位置（N 个点）
xHalf = x(1:N-1) + 0.5 * dx; % v 的位置（N-1 个点）

pulse = @(s) amp * exp(-(s.^2) / (2 * sigma^2));

[pTaylor, ~] = propagate(taylorCoeff, x, xHalf, x0, sigma, amp, c, dt, nt);
[pOpt,    ~] = propagate(optCoeff,    x, xHalf, x0, sigma, amp, c, dt, nt);

pExact = pulse(x - x0 - c * elapsed);

errTaylor = pTaylor - pExact;
errOpt    = pOpt    - pExact;

windowMask = x >= viewWin(1) & x <= viewWin(2);

l2Taylor   = sqrt(sum(errTaylor(windowMask).^2) / sum(windowMask));
l2Opt      = sqrt(sum(errOpt(windowMask).^2)    / sum(windowMask));
maxTaylor  = max(abs(errTaylor(windowMask)));
maxOpt     = max(abs(errOpt(windowMask)));

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1280, 480]);
tl = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

plot_panel(nexttile, x, pTaylor, pExact, viewWin, ...
    '(a)', 'Taylor FD', 'Exact solution');
plot_panel(nexttile, x, pOpt, pExact, viewWin, ...
    '(b)', 'Optimized FD', 'Exact solution');

outFile = fullfile(projectDir, 'results_standard', 'figures', 'exp_1d_pulse_vs_exact.png');
ewm_save_figure(fig, outFile);
close(fig);

fprintf('已保存：%s\n', outFile);
fprintf('网格 N = %d，dx = %.3f，c = %.3f，dt = %.5g（CFL = %.3f）\n', ...
    N, dx, c, dt, c * dt / dx);
fprintf('传播步数 nt = %d，传播时间 = %.4f；脉冲从 x/Δx = %g 移到 %g\n', ...
    nt, elapsed, x0, x0 + c * elapsed / dx);
fprintf('Taylor 系数：[%s]\n', strjoin(arrayfun(@(v) sprintf('%.6f', v), taylorCoeff, 'UniformOutput', false), ', '));
fprintf('优化系数 ：[%s]\n', strjoin(arrayfun(@(v) sprintf('%.6f', v), optCoeff,    'UniformOutput', false), ', '));
fprintf('窗口 x/Δx ∈ [%g, %g] 内误差：\n', viewWin(1), viewWin(2));
fprintf('  Taylor   ：L2 = %.4g，最大|误差| = %.4g\n', l2Taylor, maxTaylor);
fprintf('  优化     ：L2 = %.4g，最大|误差| = %.4g\n', l2Opt,    maxOpt);
fprintf('  优化/Taylor 比值：L2 = %.3f，最大 = %.3f\n', l2Opt/l2Taylor, maxOpt/maxTaylor);

% ============== 子函数 ==============

function [p, v] = propagate(coeff, x, xHalf, x0, sigma, amp, c, dt, nt)
% 交错网格 1D 声学方程：dp/dt = -dv/dx，dv/dt = -dp/dx（ρ = c = 1）
% p 在整数节点 x（N 个），v 在半节点 xHalf（N-1 个，位于 x(i) 与 x(i+1) 之间）。
order = numel(coeff);
N = numel(x);
NH = numel(xHalf);  % = N - 1

pulse = @(s) amp * exp(-(s.^2) / (2 * sigma^2));

p = pulse(x - x0).';
v = pulse(xHalf - x0 + c * dt / 2).';   % v 初始时刻 t = -dt/2，对应实参 +c·dt/2

for it = 1:nt
    dpdx = zeros(NH, 1);
    for n = 1:order
        idx = n : (NH - n + 1);          % v 索引：保证 p(idx+n)、p(idx-n+1) 合法
        dpdx(idx) = dpdx(idx) + coeff(n) * (p(idx + n) - p(idx - n + 1));
    end
    v = v - dt * dpdx;

    dvdx = zeros(N, 1);
    for n = 1:order
        idx = (n + 1) : (N - n);         % p 索引：保证 v(idx+n-1)、v(idx-n) 合法
        dvdx(idx) = dvdx(idx) + coeff(n) * (v(idx + n - 1) - v(idx - n));
    end
    p = p - dt * dvdx;
end

p = p.';
v = v.';
end

function plot_panel(ax, x, pFD, pExact, viewWin, tag, fdLabel, exactLabel)
hold(ax, 'on');
plot(ax, x, pFD,    'k-',  'LineWidth', 1.3, 'DisplayName', fdLabel);
plot(ax, x, pExact, 'k--', 'LineWidth', 1.0, 'DisplayName', exactLabel);
hold(ax, 'off');
grid(ax, 'on');
box(ax, 'on');
xTickStart = viewWin(1) - mod(viewWin(1), 5);
set(ax, 'XLim', viewWin, 'YLim', [-0.1, 0.6], ...
    'XTick', xTickStart:5:viewWin(2), 'Layer', 'top', 'LineWidth', 0.8);
xlabel(ax, 'x / \Deltax');
ylabel(ax, 'Amplitude');
legend(ax, 'Location', 'northwest', 'Box', 'on');
text(ax, viewWin(1) + 1, 0.55, tag, 'FontWeight', 'bold', 'FontSize', 12);
end
