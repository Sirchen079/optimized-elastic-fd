function [ratio, boundaryEnergy, totalEnergy] = ewm_boundary_ratio(vx, vz, width)
%EWM_BOUNDARY_RATIO 计算外侧边框区域的能量占比。

energy = vx .^ 2 + vz .^ 2;
totalEnergy = sum(energy(:));

mask = false(size(energy));
mask(1:width, :) = true;
mask(end-width+1:end, :) = true;
mask(:, 1:width) = true;
mask(:, end-width+1:end) = true;

boundaryEnergy = sum(energy(mask));
ratio = boundaryEnergy / (totalEnergy + eps);
end
