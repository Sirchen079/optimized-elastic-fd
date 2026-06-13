function [ddx, ddz] = ewm_pml2d(vp, nLayer, dx, dz)
%EWM_PML2D 基于可信 pml2d.m 的 PML 衰减剖面。

[NZ, NX] = size(vp);
ddz = zeros(NZ, 1);
ddx = zeros(1, NX);

if nLayer <= 0
    return;
end

vpMax = max(vp(:));
Lx = nLayer * dx;
Lz = nLayer * dz;
R = 0.001 * 10 ^ (-(log10(nLayer) - 1) / log10(2) - 3);
d0x = -log(R) * (3 * vpMax) / (2 * Lx);
d0z = -log(R) * (3 * vpMax) / (2 * Lz);

for i = 1:nLayer
    wx = ((nLayer - i) / nLayer) ^ 2;
    wz = wx;
    ddz(i) = d0z * wz;
    ddz(NZ - i + 1) = ddz(i);
    ddx(i) = d0x * wx;
    ddx(NX - i + 1) = ddx(i);
end
end
