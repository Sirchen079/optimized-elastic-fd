function cmap = ewm_wavefield_colormap(n)
%EWM_WAVEFIELD_COLORMAP 蓝-白-红发散色图，适合波场正负振幅可视化。
if nargin < 1
    n = 256;
end
half = ceil(n / 2);
r = [linspace(0.10, 1.00, half), linspace(1.00, 0.80, n - half)]';
g = [linspace(0.20, 1.00, half), linspace(1.00, 0.15, n - half)]';
b = [linspace(0.72, 1.00, half), linspace(1.00, 0.10, n - half)]';
cmap = [r, g, b];
end
