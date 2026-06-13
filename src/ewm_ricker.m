function src = ewm_ricker(t, f0, delayCycles)
%EWM_RICKER 主频为 f0 的 Ricker 子波。

if nargin < 3 || isempty(delayCycles)
    delayCycles = 1.5;
end

t0 = delayCycles / f0;
a = pi * f0 * (t - t0);
src = (1 - 2 * a.^2) .* exp(-a.^2);
end
