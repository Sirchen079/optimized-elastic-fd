%PLOT_DISPERSION_COMPARISON 独立脚本：画师兄风格的交错网格一阶导数频散误差图。
%
% 用法：在 MATLAB 里直接运行本文件，无需先跑完整仿真。
% 脚本自动从 results_standard/coefficients.txt 读取已保存的优化系数，
% 同时计算同阶 Taylor 系数作为对比，最终输出 dispersion_signed.png。
%
% 如果 coefficients.txt 不存在，脚本会跳过文件读取，用默认硬编码系数作为演示。

clear; clc;
projectDir = fileparts(mfilename('fullpath'));
addpath(fullfile(projectDir, 'src'));

% -------- 1. 读取优化系数 --------
coeffFile = fullfile(projectDir, 'results_standard', 'coefficients.txt');
optimizedCoeff = [];
if isfile(coeffFile)
    fid = fopen(coeffFile, 'r');
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if startsWith(line, '优化系数 =')
            nums = sscanf(line(length('优化系数 =') + 1:end), '%f');
            optimizedCoeff = nums(:).';
            break;
        end
    end
    fclose(fid);
end

if isempty(optimizedCoeff)
    warning('未找到 coefficients.txt，使用内置默认 10 阶优化系数。');
    optimizedCoeff = [1.23717720156044, -0.108727049806825, ...
                      0.0237950229290865, -0.0052569485184926, 0.000758608093783266];
end

order = numel(optimizedCoeff);   % 单侧系数个数，对应 2*order 阶

% -------- 2. 计算同阶 Taylor 系数 --------
standardCoeff = ewm_fd_coefficients('standard', order);

fprintf('Taylor   系数 (%d阶): %s\n', 2*order, sprintf('%.6g  ', standardCoeff));
fprintf('优化系数 (%d阶): %s\n', 2*order, sprintf('%.6g  ', optimizedCoeff));

% -------- 3. 从配置读取目标带宽和误差阈值 --------
khMax       = 0.60 * pi;   % 与 ewm_default_config 一致
targetError = 1e-4;

% -------- 4. 画图并保存 --------
outFile = fullfile(projectDir, 'dispersion_signed.png');
ewm_plot_dispersion_signed(standardCoeff, optimizedCoeff, khMax, targetError, outFile);
fprintf('\n图已保存到：%s\n', outFile);
imshow(outFile);   % 直接在 MATLAB 图窗里预览
