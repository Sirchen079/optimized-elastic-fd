function xl = ewm_wave_xlim(xKm, varargin)
%EWM_WAVE_XLIM 计算波场快照的水平有效显示范围，用于裁掉左右空白区域。
%
% 输入：
%   xKm       水平坐标向量（单位 km），长度为 nx。
%   varargin  一个或多个波场数组，尺寸为 nz x nx 或 nz x nx x nSnap。
%
% 输出：
%   xl        [xmin, xmax]，覆盖所有传入波场中振幅显著的列，并留少量边距。
%
% 原理：对每一列求跨深度（及跨快照）的能量，归一化后取超过阈值的列范围，
%       两端各留出约 6%% 跨度的边距，避免把波前裁到面板边缘。

energy = [];
for i = 1:numel(varargin)
    f = varargin{i};
    f = reshape(f, size(f, 1), size(f, 2), []);
    e = squeeze(sqrt(sum(sum(f.^2, 1), 3)));   % nx x 1 的逐列能量
    e = e(:);
    if isempty(energy)
        energy = e;
    else
        energy = energy + e;
    end
end

nx = numel(xKm);
if isempty(energy) || max(energy) <= 0
    xl = [xKm(1), xKm(end)];
    return;
end

energy = energy / max(energy);
idx = find(energy > 0.01);
if isempty(idx)
    xl = [xKm(1), xKm(end)];
    return;
end

i1 = idx(1);
i2 = idx(end);
span = i2 - i1;
pad = max(round(0.06 * span), round(0.02 * nx));
i1 = max(1, i1 - pad);
i2 = min(nx, i2 + pad);
xl = [xKm(i1), xKm(i2)];
end
