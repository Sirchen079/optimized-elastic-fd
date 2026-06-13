function result = ewm_staggered_solver(model, sim, coeff, usePml, tag)
%EWM_STAGGERED_SOLVER 二维弹性交错网格有限差分求解器。
%
% 导数模板采用标准二维弹性交错网格有限差分格式（含 PML 阻尼），
% 差分系数通过参数传入。

if nargin < 5 || isempty(tag)
    tag = 'staggered';
end

coeff = coeff(:).';
M = numel(coeff);
if usePml
    nLayer = sim.nPml;
else
    % 仅使用非吸收边界的虚拟层以支撑高阶差分模板，
    % 输出时仍裁剪回原始物理模型区域。
    nLayer = M;
end

vp = ewm_pad2d(model.vp, nLayer);
vs = ewm_pad2d(model.vs, nLayer);
rho = ewm_pad2d(model.rho, nLayer);
[NZ, NX] = size(vp);

ep = ewm_prepare_elastic_params(vp, vs, rho);

iz = (M+1):(NZ-M);
ix = (M+1):(NX-M);
dt = sim.dt;
dx = model.dx;
dz = model.dz;

t = (0:sim.nt-1) * dt;
src = sim.sourceAmplitude * ewm_ricker(t, sim.f0, sim.sourceDelayCycles);

% 震源类型：
%   'explosive'   (默认) 各向同性应力源，加在 τ_xx、τ_zz 上，只辐射 P
%   'pointForceX' 水平点力，加在 vx 上；辐射 P+SV，可用于演示 S 频散
%   'pointForceZ' 垂直点力，加在 vz 上
if isfield(sim, 'sourceType') && ~isempty(sim.sourceType)
    sourceType = lower(sim.sourceType);
else
    sourceType = 'explosive';
end

srcI0 = min(max(round(sim.sourceDepthM / dz) + 1, M+1), model.nz - M);
srcJ0 = min(max(round(1 + (model.nx - 1) * sim.sourceXFraction), M+1), model.nx - M);
srcI = srcI0 + nLayer;
srcJ = srcJ0 + nLayer;

% 震源空间分布：默认单点（stamp = 1），可选高斯加宽以压制 δ 源在网格上
% 注入的高 kh 噪声。stamp 应为奇数边长的方阵，且 sum(stamp(:)) == 1。
if isfield(sim, 'sourceStamp') && ~isempty(sim.sourceStamp)
    stamp = sim.sourceStamp;
else
    stamp = 1;
end
stampHalf = (size(stamp, 1) - 1) / 2;
stampI = srcI + (-stampHalf:stampHalf);
stampJ = srcJ + (-stampHalf:stampHalf);

vx = zeros(NZ, NX);
vz = zeros(NZ, NX);
tau_xx = zeros(NZ, NX);
tau_zz = zeros(NZ, NX);
tau_xz = zeros(NZ, NX);

if usePml
    vx_x = zeros(NZ, NX);
    vx_z = zeros(NZ, NX);
    vz_x = zeros(NZ, NX);
    vz_z = zeros(NZ, NX);
    tau_xx_x = zeros(NZ, NX);
    tau_xx_z = zeros(NZ, NX);
    tau_zz_x = zeros(NZ, NX);
    tau_zz_z = zeros(NZ, NX);
    tau_xz_x = zeros(NZ, NX);
    tau_xz_z = zeros(NZ, NX);

    [ddx, ddz] = ewm_pml2d(vp, nLayer, dx, dz);
    coefx1 = (1 - ddx * dt / 2) ./ (1 + ddx * dt / 2);
    coefx2 = dt ./ (1 + ddx * dt / 2);
    coefz1 = (1 - ddz * dt / 2) ./ (1 + ddz * dt / 2);
    coefz2 = dt ./ (1 + ddz * dt / 2);
end

snapshotIdx = ewm_pick_snapshot_indices(sim, dt);
snapshots = zeros(model.nz, model.nx, numel(snapshotIdx));
snapTimes = t(snapshotIdx);
snapCounter = 0;

% 可选：在指定接收点逐时间步记录 vz、vx 单道时序（用于波形对比）。
recordTrace = isfield(sim, 'receivers') && ~isempty(sim.receivers);
if recordTrace
    recIz = sim.receivers.iz(:);
    recJx = sim.receivers.ix(:);
    nRec = numel(recIz);
    traceVz = zeros(sim.nt, nRec);
    traceVx = zeros(sim.nt, nRec);
end

energyStep = max(1, sim.energyStride);
energyCount = ceil(sim.nt / energyStep);
energyTime = zeros(energyCount, 1);
boundaryRatio = zeros(energyCount, 1);
totalVelocityEnergy = zeros(energyCount, 1);
boundaryVelocityEnergy = zeros(energyCount, 1);
energyCursor = 0;
edgeWidth = max(M + 3, round(0.08 * min(model.nz, model.nx)));

fprintf('  正在运行：%s\n', ewm_solver_display_name(tag));
progressEvery = max(1, floor(sim.nt / 5));

for it = 1:sim.nt
    if usePml
        dtau_xxdx = sg_dx_forward(tau_xx, iz, ix, coeff, dx);
        dtau_zzdz = sg_dz_forward(tau_zz, iz, ix, coeff, dz);
        dtau_xzdx = sg_dx_backward(tau_xz, iz, ix, coeff, dx);
        dtau_xzdz = sg_dz_backward(tau_xz, iz, ix, coeff, dz);

        vx_x(iz, ix) = vx_x(iz, ix) .* coefx1(ix) + ...
            ep.invrhox(iz, ix) .* dtau_xxdx .* coefx2(ix);
        vx_z(iz, ix) = vx_z(iz, ix) .* coefz1(iz) + ...
            ep.invrhox(iz, ix) .* dtau_xzdz .* coefz2(iz);
        vz_x(iz, ix) = vz_x(iz, ix) .* coefx1(ix) + ...
            ep.invrhoz(iz, ix) .* dtau_xzdx .* coefx2(ix);
        vz_z(iz, ix) = vz_z(iz, ix) .* coefz1(iz) + ...
            ep.invrhoz(iz, ix) .* dtau_zzdz .* coefz2(iz);
        vx = vx_x + vx_z;
        vz = vz_x + vz_z;

        dvxdx = sg_dx_backward(vx, iz, ix, coeff, dx);
        dvxdz = sg_dz_forward(vx, iz, ix, coeff, dz);
        dvzdx = sg_dx_forward(vz, iz, ix, coeff, dx);
        dvzdz = sg_dz_backward(vz, iz, ix, coeff, dz);

        tau_xx_x(iz, ix) = tau_xx_x(iz, ix) .* coefx1(ix) + ...
            ep.c11(iz, ix) .* dvxdx .* coefx2(ix);
        tau_xx_z(iz, ix) = tau_xx_z(iz, ix) .* coefz1(iz) + ...
            ep.c13(iz, ix) .* dvzdz .* coefz2(iz);
        tau_zz_x(iz, ix) = tau_zz_x(iz, ix) .* coefx1(ix) + ...
            ep.c13(iz, ix) .* dvxdx .* coefx2(ix);
        tau_zz_z(iz, ix) = tau_zz_z(iz, ix) .* coefz1(iz) + ...
            ep.c33(iz, ix) .* dvzdz .* coefz2(iz);
        tau_xz_x(iz, ix) = tau_xz_x(iz, ix) .* coefx1(ix) + ...
            ep.c44(iz, ix) .* dvzdx .* coefx2(ix);
        tau_xz_z(iz, ix) = tau_xz_z(iz, ix) .* coefz1(iz) + ...
            ep.c44(iz, ix) .* dvxdz .* coefz2(iz);

        switch sourceType
            case 'pointforcex'
                inj = 0.5 * dt * src(it) * stamp;
                vx_x(stampI, stampJ) = vx_x(stampI, stampJ) + inj;
                vx_z(stampI, stampJ) = vx_z(stampI, stampJ) + inj;
                vx = vx_x + vx_z;
            case 'pointforcez'
                inj = 0.5 * dt * src(it) * stamp;
                vz_x(stampI, stampJ) = vz_x(stampI, stampJ) + inj;
                vz_z(stampI, stampJ) = vz_z(stampI, stampJ) + inj;
                vz = vz_x + vz_z;
            otherwise   % 'explosive'
                halfAmp = 0.5 * src(it) * stamp;
                tau_xx_x(stampI, stampJ) = tau_xx_x(stampI, stampJ) + halfAmp;
                tau_xx_z(stampI, stampJ) = tau_xx_z(stampI, stampJ) + halfAmp;
                tau_zz_x(stampI, stampJ) = tau_zz_x(stampI, stampJ) + halfAmp;
                tau_zz_z(stampI, stampJ) = tau_zz_z(stampI, stampJ) + halfAmp;
        end

        tau_xx = tau_xx_x + tau_xx_z;
        tau_zz = tau_zz_x + tau_zz_z;
        tau_xz = tau_xz_x + tau_xz_z;
    else
        dtau_xxdx = sg_dx_forward(tau_xx, iz, ix, coeff, dx);
        dtau_zzdz = sg_dz_forward(tau_zz, iz, ix, coeff, dz);
        dtau_xzdx = sg_dx_backward(tau_xz, iz, ix, coeff, dx);
        dtau_xzdz = sg_dz_backward(tau_xz, iz, ix, coeff, dz);

        vx(iz, ix) = vx(iz, ix) + dt * ep.invrhox(iz, ix) .* ...
            (dtau_xxdx + dtau_xzdz);
        vz(iz, ix) = vz(iz, ix) + dt * ep.invrhoz(iz, ix) .* ...
            (dtau_xzdx + dtau_zzdz);

        dvxdx = sg_dx_backward(vx, iz, ix, coeff, dx);
        dvxdz = sg_dz_forward(vx, iz, ix, coeff, dz);
        dvzdx = sg_dx_forward(vz, iz, ix, coeff, dx);
        dvzdz = sg_dz_backward(vz, iz, ix, coeff, dz);

        tau_xx(iz, ix) = tau_xx(iz, ix) + dt * ( ...
            ep.c11(iz, ix) .* dvxdx + ep.c13(iz, ix) .* dvzdz);
        tau_zz(iz, ix) = tau_zz(iz, ix) + dt * ( ...
            ep.c13(iz, ix) .* dvxdx + ep.c33(iz, ix) .* dvzdz);
        tau_xz(iz, ix) = tau_xz(iz, ix) + dt * ...
            ep.c44(iz, ix) .* (dvzdx + dvxdz);

        switch sourceType
            case 'pointforcex'
                vx(stampI, stampJ) = vx(stampI, stampJ) + dt * src(it) * stamp;
            case 'pointforcez'
                vz(stampI, stampJ) = vz(stampI, stampJ) + dt * src(it) * stamp;
            otherwise   % 'explosive'
                tau_xx(stampI, stampJ) = tau_xx(stampI, stampJ) + src(it) * stamp;
                tau_zz(stampI, stampJ) = tau_zz(stampI, stampJ) + src(it) * stamp;
        end
    end

    if any(snapshotIdx == it)
        snapCounter = snapCounter + 1;
        snapshots(:, :, snapCounter) = vz(nLayer+1:nLayer+model.nz, ...
            nLayer+1:nLayer+model.nx);
    end

    if recordTrace
        for r = 1:nRec
            traceVz(it, r) = vz(recIz(r) + nLayer, recJx(r) + nLayer);
            traceVx(it, r) = vx(recIz(r) + nLayer, recJx(r) + nLayer);
        end
    end

    if mod(it - 1, energyStep) == 0
        energyCursor = energyCursor + 1;
        vxPhys = vx(nLayer+1:nLayer+model.nz, nLayer+1:nLayer+model.nx);
        vzPhys = vz(nLayer+1:nLayer+model.nz, nLayer+1:nLayer+model.nx);
        energyTime(energyCursor) = t(it);
        [boundaryRatio(energyCursor), boundaryVelocityEnergy(energyCursor), ...
            totalVelocityEnergy(energyCursor)] = ewm_boundary_ratio(vxPhys, vzPhys, edgeWidth);
    end

    if mod(it, progressEvery) == 0 || it == sim.nt
        fprintf('    时间步 %d / %d\n', it, sim.nt);
    end
end

result = struct();
result.tag = tag;
result.scheme = 'staggered';
result.usePml = usePml;
result.coeff = coeff;
result.dt = dt;
result.nt = sim.nt;
result.snapshots.vz = snapshots;
result.snapshots.time = snapTimes;
result.energy.time = energyTime(1:energyCursor);
result.energy.boundaryRatio = boundaryRatio(1:energyCursor);
result.energy.totalVelocity = totalVelocityEnergy(1:energyCursor);
result.energy.boundaryVelocity = boundaryVelocityEnergy(1:energyCursor);
result.source.index = [srcI0, srcJ0];
result.source.depthM = (srcI0 - 1) * dz;
result.source.xM = (srcJ0 - 1) * dx;

if recordTrace
    result.trace.time = t(:);
    result.trace.vz = traceVz;
    result.trace.vx = traceVx;
    result.trace.iz = recIz;
    result.trace.ix = recJx;
    result.trace.zM = (recIz - 1) * dz;
    result.trace.xM = (recJx - 1) * dx;
end

if isfield(sim, 'outputDir') && ~isempty(sim.outputDir)
    save(fullfile(sim.outputDir, [tag, '.mat']), 'result');
end
end

function d = sg_dx_forward(a, iz, ix, coeff, dx)
d = zeros(numel(iz), numel(ix));
for m = 1:numel(coeff)
    d = d + coeff(m) * (a(iz, ix + m) - a(iz, ix - m + 1));
end
d = d / dx;
end

function d = sg_dx_backward(a, iz, ix, coeff, dx)
d = zeros(numel(iz), numel(ix));
for m = 1:numel(coeff)
    d = d + coeff(m) * (a(iz, ix + m - 1) - a(iz, ix - m));
end
d = d / dx;
end

function d = sg_dz_forward(a, iz, ix, coeff, dz)
d = zeros(numel(iz), numel(ix));
for m = 1:numel(coeff)
    d = d + coeff(m) * (a(iz + m, ix) - a(iz - m + 1, ix));
end
d = d / dz;
end

function d = sg_dz_backward(a, iz, ix, coeff, dz)
d = zeros(numel(iz), numel(ix));
for m = 1:numel(coeff)
    d = d + coeff(m) * (a(iz + m - 1, ix) - a(iz - m, ix));
end
d = d / dz;
end
