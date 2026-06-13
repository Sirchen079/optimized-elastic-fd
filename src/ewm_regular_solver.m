function result = ewm_regular_solver(model, sim, tag)
%EWM_REGULAR_SOLVER 常规同位网格弹性有限差分基线求解器。
%
% 所有速度和应力分量均存储在同一网格节点上，
% 空间导数采用二阶中心差分。

if nargin < 3 || isempty(tag)
    tag = 'regular';
end

vp = model.vp;
vs = model.vs;
rho = model.rho;
[nz, nx] = size(vp);

mu = rho .* vs .^ 2;
lambda = rho .* (vp .^ 2 - 2 * vs .^ 2);
c11 = lambda + 2 * mu;
c33 = c11;
c13 = lambda;
c44 = mu;

dt = sim.dt;
dx = model.dx;
dz = model.dz;

t = (0:sim.nt-1) * dt;
src = sim.sourceAmplitude * ewm_ricker(t, sim.f0, sim.sourceDelayCycles);

srcI = min(max(round(sim.sourceDepthM / dz) + 1, 2), nz - 1);
srcJ = min(max(round(1 + (nx - 1) * sim.sourceXFraction), 2), nx - 1);

vx = zeros(nz, nx);
vz = zeros(nz, nx);
tau_xx = zeros(nz, nx);
tau_zz = zeros(nz, nx);
tau_xz = zeros(nz, nx);

iz = 2:nz-1;
ix = 2:nx-1;

snapshotIdx = ewm_pick_snapshot_indices(sim, dt);
snapshots = zeros(nz, nx, numel(snapshotIdx));
snapTimes = t(snapshotIdx);
snapCounter = 0;

energyStep = max(1, sim.energyStride);
energyCount = ceil(sim.nt / energyStep);
energyTime = zeros(energyCount, 1);
boundaryRatio = zeros(energyCount, 1);
totalVelocityEnergy = zeros(energyCount, 1);
boundaryVelocityEnergy = zeros(energyCount, 1);
energyCursor = 0;
edgeWidth = max(8, round(0.08 * min(nz, nx)));

fprintf('  正在运行：%s\n', ewm_solver_display_name(tag));
progressEvery = max(1, floor(sim.nt / 5));

for it = 1:sim.nt
    dtau_xxdx = ddx_center(tau_xx, iz, ix, dx);
    dtau_zzdz = ddz_center(tau_zz, iz, ix, dz);
    dtau_xzdx = ddx_center(tau_xz, iz, ix, dx);
    dtau_xzdz = ddz_center(tau_xz, iz, ix, dz);

    vx(iz, ix) = vx(iz, ix) + dt ./ rho(iz, ix) .* (dtau_xxdx + dtau_xzdz);
    vz(iz, ix) = vz(iz, ix) + dt ./ rho(iz, ix) .* (dtau_xzdx + dtau_zzdz);

    dvxdx = ddx_center(vx, iz, ix, dx);
    dvxdz = ddz_center(vx, iz, ix, dz);
    dvzdx = ddx_center(vz, iz, ix, dx);
    dvzdz = ddz_center(vz, iz, ix, dz);

    tau_xx(iz, ix) = tau_xx(iz, ix) + dt * ( ...
        c11(iz, ix) .* dvxdx + c13(iz, ix) .* dvzdz);
    tau_zz(iz, ix) = tau_zz(iz, ix) + dt * ( ...
        c13(iz, ix) .* dvxdx + c33(iz, ix) .* dvzdz);
    tau_xz(iz, ix) = tau_xz(iz, ix) + dt * ...
        c44(iz, ix) .* (dvzdx + dvxdz);

    tau_xx(srcI, srcJ) = tau_xx(srcI, srcJ) + src(it);
    tau_zz(srcI, srcJ) = tau_zz(srcI, srcJ) + src(it);

    if any(snapshotIdx == it)
        snapCounter = snapCounter + 1;
        snapshots(:, :, snapCounter) = vz;
    end

    if mod(it - 1, energyStep) == 0
        energyCursor = energyCursor + 1;
        energyTime(energyCursor) = t(it);
        [boundaryRatio(energyCursor), boundaryVelocityEnergy(energyCursor), ...
            totalVelocityEnergy(energyCursor)] = ewm_boundary_ratio(vx, vz, edgeWidth);
    end

    if mod(it, progressEvery) == 0 || it == sim.nt
        fprintf('    时间步 %d / %d\n', it, sim.nt);
    end
end

result = struct();
result.tag = tag;
result.scheme = 'regular';
result.usePml = false;
result.dt = dt;
result.nt = sim.nt;
result.snapshots.vz = snapshots;
result.snapshots.time = snapTimes;
result.energy.time = energyTime(1:energyCursor);
result.energy.boundaryRatio = boundaryRatio(1:energyCursor);
result.energy.totalVelocity = totalVelocityEnergy(1:energyCursor);
result.energy.boundaryVelocity = boundaryVelocityEnergy(1:energyCursor);
result.source.index = [srcI, srcJ];
result.source.depthM = (srcI - 1) * dz;
result.source.xM = (srcJ - 1) * dx;

if isfield(sim, 'outputDir') && ~isempty(sim.outputDir)
    save(fullfile(sim.outputDir, [tag, '.mat']), 'result');
end
end

function d = ddx_center(a, iz, ix, dx)
d = (a(iz, ix + 1) - a(iz, ix - 1)) / (2 * dx);
end

function d = ddz_center(a, iz, ix, dz)
d = (a(iz + 1, ix) - a(iz - 1, ix)) / (2 * dz);
end
