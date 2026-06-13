function padded = ewm_pad2d(modelParam, nLayer)
%EWM_PAD2D 通过复制边缘值扩展二维模型。

if nLayer <= 0
    padded = modelParam;
    return;
end

[nz, nx] = size(modelParam);
NZ = nz + 2 * nLayer;
NX = nx + 2 * nLayer;

padded = zeros(NZ, NX);
padded(nLayer+1:nLayer+nz, nLayer+1:nLayer+nx) = modelParam;
padded(1:nLayer, :) = padded(nLayer+1, :) .* ones(nLayer, NX);
padded(nLayer+nz+1:NZ, :) = padded(nLayer+nz, :) .* ones(nLayer, NX);
padded(:, 1:nLayer) = padded(:, nLayer+1) .* ones(NZ, nLayer);
padded(:, nLayer+nx+1:NX) = padded(:, nLayer+nx) .* ones(NZ, nLayer);
end
