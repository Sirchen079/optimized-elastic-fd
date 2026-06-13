function dt = ewm_stable_dt(model, coeff, cfl)
%EWM_STABLE_DT 显式弹性求解器的保守稳定时间步长。

if nargin < 3 || isempty(cfl)
    cfl = 0.42;
end

maxVp = max(model.vp(:));
dxMin = min(model.dx, model.dz);
stencilGain = max(1, sum(abs(coeff)));
dt = cfl * dxMin / (sqrt(2) * maxVp * stencilGain);
end
