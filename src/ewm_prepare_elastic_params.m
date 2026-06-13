function ep = ewm_prepare_elastic_params(vp, vs, rho)
%EWM_PREPARE_ELASTIC_PARAMS 构建交错网格弹性参数。

[NZ, NX] = size(vp);

mu = rho .* vs .^ 2;
lambda = rho .* (vp .^ 2 - 2 * vs .^ 2);

ep = struct();
ep.c11 = lambda + 2 * mu;
ep.c13 = lambda;
ep.c33 = ep.c11;

ep.invrhox = 1 ./ rho;
ep.invrhoz = 1 ./ rho;
ep.invrhox(:, 1:NX-1) = 2 ./ (rho(:, 1:NX-1) + rho(:, 2:NX));
ep.invrhox(:, NX) = ep.invrhox(:, NX-1);
ep.invrhoz(1:NZ-1, :) = 2 ./ (rho(1:NZ-1, :) + rho(2:NZ, :));
ep.invrhoz(NZ, :) = ep.invrhoz(NZ-1, :);

ep.c44 = mu;
ep.c44(1:NZ-1, 1:NX-1) = 0.25 * ( ...
    mu(1:NZ-1, 1:NX-1) + mu(1:NZ-1, 2:NX) + ...
    mu(2:NZ, 1:NX-1) + mu(2:NZ, 2:NX));
ep.c44(NZ, :) = ep.c44(NZ-1, :);
ep.c44(:, NX) = ep.c44(:, NX-1);
end
