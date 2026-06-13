function [modelRef, cropInfo, simRef] = ewm_make_boundary_reference_case(model, sim, extraGrid)
%EWM_MAKE_BOUNDARY_REFERENCE_CASE 扩展计算域使边界远离震源区域。

if nargin < 3 || isempty(extraGrid)
    extraGrid = 25;
end

modelRef = model;
modelRef.vp = ewm_pad2d(model.vp, extraGrid);
modelRef.vs = ewm_pad2d(model.vs, extraGrid);
modelRef.rho = ewm_pad2d(model.rho, extraGrid);
[modelRef.nz, modelRef.nx] = size(modelRef.vp);
modelRef.x = (0:modelRef.nx-1) * modelRef.dx;
modelRef.z = (0:modelRef.nz-1) * modelRef.dz;
modelRef.source = sprintf('%s, extended by %d grid points', model.source, extraGrid);

cropInfo = struct();
cropInfo.z = (extraGrid + 1):(extraGrid + model.nz);
cropInfo.x = (extraGrid + 1):(extraGrid + model.nx);
cropInfo.extraGrid = extraGrid;

sourceXIndex = round(1 + (model.nx - 1) * sim.sourceXFraction);
sourceXIndexRef = extraGrid + sourceXIndex;

simRef = sim;
simRef.sourceDepthM = sim.sourceDepthM + extraGrid * model.dz;
simRef.sourceXFraction = (sourceXIndexRef - 1) / (modelRef.nx - 1);
simRef.nPml = sim.nPml;
end
