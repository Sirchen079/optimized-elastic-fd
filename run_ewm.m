function results = run_ewm(mode)
%RUN_EWM 弹性波论文实验主入口。
%
% 用法：
%   run_ewm              % 快速验证（80 m 缓存）
%   run_ewm('preview')     % 同上
%   run_ewm('standard')    % 论文模式（10 m 网格，12 Hz 震源）
%   run_ewm('highfreq') % 论文高频模式（10 m 网格，25 Hz 震源）

if nargin < 1 || isempty(mode)
    mode = 'preview';
end

projectDir = fileparts(mfilename('fullpath'));
addpath(fullfile(projectDir, 'src'));
ewm_apply_chinese_style();

cfg = ewm_default_config(projectDir, mode);
ewm_ensure_dir(cfg.output.dir);
cfg.output.figuresDir = fullfile(cfg.output.dir, 'figures');
ewm_ensure_dir(cfg.output.figuresDir);
startTime = tic;
if (cfg.reference.enabled && cfg.reference.spacing == 10) || cfg.model.spacing == 10
    ewm_build_marmousi_10m_cache(cfg.model.root, false);
end

fprintf('ewm 模式：%s\n', cfg.mode);
fprintf('正在加载 Marmousi 模型：%s\n', cfg.model.root);
model = ewm_load_marmousi(cfg.model);
fprintf('模型尺寸：nz=%d, nx=%d, dz=%.3f m, dx=%.3f m\n', ...
    model.nz, model.nx, model.dz, model.dx);

fprintf('正在使用模拟退火法优化基于最大范数目标函数的差分系数...\n');
standardCoeff = ewm_fd_coefficients('standard', cfg.coeff.order);
[optimizedCoeff, optInfo] = ewm_optimize_minimax_coeffs( ...
    cfg.coeff.order, cfg.coeff.khMax, cfg.coeff.samples, cfg.coeff);

cfg.sim.dtTaylorLimit = ewm_stable_dt(model, standardCoeff, cfg.sim.cfl);
cfg.sim.dtOptimizedLimit = ewm_stable_dt(model, optimizedCoeff, cfg.sim.cfl);
cfg.sim.dt = min(cfg.sim.dtTaylorLimit, cfg.sim.dtOptimizedLimit);

assert(cfg.sim.dt <= cfg.sim.dtTaylorLimit + eps, ...
    'ewm:DtExceedsTaylorLimit', ...
    '当前 dt 超过 Taylor 系数稳定时间步长。');

assert(cfg.sim.dt <= cfg.sim.dtOptimizedLimit + eps, ...
    'ewm:DtExceedsOptimizedLimit', ...
    '当前 dt 超过优化系数稳定时间步长。');

allSnapshotTimes = [ ...
    cfg.snapshots.exp1Times(:); ...
    cfg.snapshots.exp2Times(:); ...
    cfg.snapshots.exp3Times(:) ...
];
cfg.sim.maxSnapshotTime = max(allSnapshotTimes);

cfg.sim.nt = max(cfg.sim.nt, ceil(cfg.sim.maxSnapshotTime / cfg.sim.dt) + 1);
cfg.sim.totalTime = cfg.sim.dt * (cfg.sim.nt - 1);

assert(cfg.sim.totalTime + 1e-12 >= cfg.sim.maxSnapshotTime, ...
    'ewm:SimulationTimeTooShort', ...
    '总模拟时长 %.6f s 小于最大快照时刻 %.6f s。', ...
    cfg.sim.totalTime, cfg.sim.maxSnapshotTime);

fprintf('时间步长：dt = %.12g s，nt = %d，总时长 = %.6f s，最大快照时刻 = %.6f s\n', ...
    cfg.sim.dt, cfg.sim.nt, cfg.sim.totalTime, cfg.sim.maxSnapshotTime);

ewm_plot_model(model, cfg);
ewm_plot_ricker_wavelet(cfg);

dispersion = ewm_dispersion_metrics(standardCoeff, optimizedCoeff, ...
    cfg.coeff.khMax, cfg.coeff.evalSamples);
if dispersion.optimized.maxAbsError > cfg.coeff.targetError
    error('ewm:MinimaxTargetNotMet', ...
        '模拟退火优化后的最大频散误差 %.6g 未低于目标 %.6g。', ...
        dispersion.optimized.maxAbsError, cfg.coeff.targetError);
end
ewm_plot_dispersion(standardCoeff, optimizedCoeff, dispersion, cfg, optInfo);
ewm_plot_sa_convergence(optInfo, cfg.coeff.targetError, ...
    fullfile(cfg.output.figuresDir, 'sa_convergence.png'));

fprintf('\n标准 Taylor 系数：\n');
fprintf('  %.12g', standardCoeff);
fprintf('\n基于最大范数目标函数的优化系数：\n');
fprintf('  %.12g', optimizedCoeff);
fprintf('\n波数范围 kh ∈ [0, %.3fπ] 的最大相对导数误差：\n', cfg.coeff.khMax / pi);
fprintf('  标准系数  = %.6g\n', dispersion.standard.maxAbsError);
fprintf('  优化系数  = %.6g\n', dispersion.optimized.maxAbsError);
if dispersion.optimized.maxAbsError <= cfg.coeff.targetError
    fprintf('  结论：优化误差已低于 %.1e。\n\n', cfg.coeff.targetError);
else
    fprintf('  警告：优化误差未低于 %.1e，请减小 khMax 或提高差分阶数。\n\n', cfg.coeff.targetError);
end

simBase = cfg.sim;
simBase.standardCoeff = standardCoeff;
simBase.optimizedCoeff = optimizedCoeff;

fprintf('实验1：正在对比常规网格与交错网格...\n');
simExp1 = simBase;
simExp1.snapshotTimes = cfg.snapshots.exp1Times;
regularNoAbsorb = ewm_regular_solver(model, simExp1, 'exp1_regular_noabsorb');
staggeredNoAbsorbExp1 = ewm_staggered_solver(model, simExp1, ...
    standardCoeff, false, 'exp1_staggered_noabsorb_standard');
ewm_plot_comparison(model, regularNoAbsorb, staggeredNoAbsorbExp1, ...
    '常规网格', '交错网格', ...
    fullfile(cfg.output.figuresDir, 'exp1_regular_vs_staggered.png'), 2);

fprintf('实验2：正在对比无吸收边界与 PML 吸收边界...\n');
simExp2 = simBase;
simExp2.snapshotTimes = cfg.snapshots.exp2Times;
simExp2.nt = max(simExp2.nt, ceil(max(simExp2.snapshotTimes) / simExp2.dt) + 1);
staggeredNoAbsorbExp2 = ewm_staggered_solver(model, simExp2, ...
    standardCoeff, false, 'exp2_staggered_noabsorb_standard');
staggeredPmlStandardExp2 = ewm_staggered_solver(model, simExp2, ...
    standardCoeff, true, 'exp2_staggered_pml_standard');
ewm_plot_pml_snapshot_pair(model, staggeredNoAbsorbExp2, staggeredPmlStandardExp2, ...
    fullfile(cfg.output.figuresDir, 'exp2_noabsorb_vs_pml.png'));
ewm_plot_energy(staggeredNoAbsorbExp2, staggeredPmlStandardExp2, ...
    '无吸收边界', 'PML 吸收边界', ...
    fullfile(cfg.output.figuresDir, 'exp2_energy_noabsorb_vs_pml.png'));
ewm_plot_pml_boundary_energy(staggeredNoAbsorbExp2, staggeredPmlStandardExp2, ...
    fullfile(cfg.output.figuresDir, 'exp2_pml_boundary_energy_time.png'), ...
    fullfile(cfg.output.dir, 'exp2_pml_boundary_energy_time.csv'));
boundaryReferenceComparison = [];
if cfg.boundaryReference.enabled
    fprintf('实验2补充：正在计算扩大计算域参考解，用于判断边界反射抑制效果...\n');
    [modelBoundaryRef, cropInfo, simBoundaryRef] = ewm_make_boundary_reference_case( ...
        model, simExp2, cfg.boundaryReference.extraGrid);
    boundaryReferenceResult = ewm_staggered_solver(modelBoundaryRef, simBoundaryRef, ...
        standardCoeff, true, 'exp2_boundary_large_domain_reference');
    boundaryReferenceComparison = ewm_plot_boundary_reference_error(model, ...
        staggeredNoAbsorbExp2, staggeredPmlStandardExp2, modelBoundaryRef, ...
        boundaryReferenceResult, cropInfo, ...
        fullfile(cfg.output.figuresDir, 'exp2_boundary_reference_error.png'), ...
        fullfile(cfg.output.dir, 'exp2_boundary_reference_metrics.txt'));
end

fprintf('实验3：正在对比标准系数与基于最大范数目标函数的优化系数...\n');
simExp3 = simBase;
simExp3.snapshotTimes = cfg.snapshots.exp3Times;
staggeredPmlStandardExp3 = ewm_staggered_solver(model, simExp3, ...
    standardCoeff, true, 'exp3_staggered_pml_standard');
staggeredPmlOptimized = ewm_staggered_solver(model, simExp3, ...
    optimizedCoeff, true, 'exp3_staggered_pml_minimax');
ewm_plot_comparison(model, staggeredPmlStandardExp3, staggeredPmlOptimized, ...
    '标准系数', '基于最大范数目标函数的优化系数', ...
    fullfile(cfg.output.figuresDir, 'exp3_standard_vs_minimax.png'));
exp3WavefieldDifference = ewm_plot_wavefield_triptych(model, ...
    staggeredPmlStandardExp3, staggeredPmlOptimized, ...
    'Taylor 系数', '基于最大范数目标函数的优化系数', ...
    fullfile(cfg.output.figuresDir, 'exp3_taylor_optimized_difference.png'), ...
    fullfile(cfg.output.dir, 'exp3_taylor_optimized_difference_metrics.csv'));
exp3ZoomInfo = ewm_plot_wavefield_zoom(model, ...
    staggeredPmlStandardExp3, staggeredPmlOptimized, ...
    'Taylor 系数', '基于最大范数目标函数的优化系数', ...
    fullfile(cfg.output.figuresDir, 'exp3_taylor_optimized_difference_zoom.png'), ...
    fullfile(cfg.output.dir, 'exp3_taylor_optimized_difference_zoom.txt'));

referenceComparison = [];
if cfg.reference.enabled && cfg.reference.spacing < cfg.model.spacing
    fprintf('实验3补充：正在计算更细网格参考解，用于定量判断标准系数与基于最大范数目标函数的优化系数哪个更接近参考结果...\n');
    refModelCfg = cfg.model;
    refModelCfg.spacing = cfg.reference.spacing;
    modelRef = ewm_load_marmousi(refModelCfg);
    simRef = simExp3;
    simRef.dt = ewm_stable_dt(modelRef, standardCoeff, cfg.sim.cfl);
    totalTime = cfg.sim.dt * (cfg.sim.nt - 1);
    simRef.nt = ceil(totalTime / simRef.dt) + 1;
    simRef.nPml = max(simExp3.nPml, round(simExp3.nPml * model.dx / modelRef.dx));
    simRef.outputDir = cfg.output.dir;
    referenceResult = ewm_staggered_solver(modelRef, simRef, ...
        standardCoeff, true, 'exp3_reference_finer_grid');
    referenceComparison = ewm_plot_reference_error(model, ...
        staggeredPmlStandardExp3, staggeredPmlOptimized, modelRef, referenceResult, ...
        fullfile(cfg.output.figuresDir, 'exp3_reference_error_comparison.png'), ...
        fullfile(cfg.output.dir, 'exp3_reference_error_metrics.txt'));
    referenceComparison.referenceNz = modelRef.nz;
    referenceComparison.referenceNx = modelRef.nx;
    referenceComparison.referenceDt = simRef.dt;
    referenceComparison.referenceNt = simRef.nt;
    referenceComparison.referencePmlLayers = simRef.nPml;
    referenceComparison.referenceSource = modelRef.source;
    referenceComparison.referenceCacheSpacing = modelRef.cacheSpacing;
    referenceComparison.referenceMaxVp = max(modelRef.vp(:));
    if isfield(modelRef, 'sourceSpacing')
        referenceComparison.referenceSourceSpacing = modelRef.sourceSpacing;
    end
    if isfield(modelRef, 'resample')
        referenceComparison.referenceResample = modelRef.resample;
    end

    targetTimes = cfg.snapshots.exp3Times(:);
    standardActual = round(targetTimes / simExp3.dt) * simExp3.dt;
    optimizedActual = standardActual;
    referenceActual = round(targetTimes / simRef.dt) * simRef.dt;
    mismatch = max(abs([standardActual - targetTimes, referenceActual - targetTimes]), [], 2);

    T = table(targetTimes, standardActual, optimizedActual, referenceActual, mismatch, ...
        'VariableNames', {'target_time','standard_actual_time','optimized_actual_time', ...
                          'reference_actual_time','max_abs_time_mismatch'});
    writetable(T, fullfile(cfg.output.dir, 'exp3_snapshot_time_alignment.csv'));

    if max(mismatch) > 0.5 * simExp3.dt
        warning('ewm:SnapshotTimeMismatch', ...
            '参考解与粗网格快照时间差较大，请检查快照索引。');
    end
end

requiredTables = ewm_write_required_tables(model, cfg, standardCoeff, ...
    optimizedCoeff, optInfo, dispersion, boundaryReferenceComparison, referenceComparison);

results = struct();
results.config = cfg;
results.model = rmfield(model, {'vp', 'vs', 'rho'});
results.coefficients.standard = standardCoeff;
results.coefficients.optimized = optimizedCoeff;
results.optimization = optInfo;
results.dispersion = dispersion;
results.regularNoAbsorb = ewm_light_result(regularNoAbsorb);
results.staggeredNoAbsorb = ewm_light_result(staggeredNoAbsorbExp2);
results.staggeredPmlStandard = ewm_light_result(staggeredPmlStandardExp2);
results.boundaryReferenceComparison = boundaryReferenceComparison;
results.staggeredPmlStandardCoeffExp3 = ewm_light_result(staggeredPmlStandardExp3);
results.staggeredPmlOptimized = ewm_light_result(staggeredPmlOptimized);
results.referenceComparison = referenceComparison;
results.exp3WavefieldDifference = exp3WavefieldDifference;
results.exp3ZoomInfo = exp3ZoomInfo;
results.requiredTables = requiredTables;

summaryFile = fullfile(cfg.output.dir, 'summary.mat');
save(summaryFile, 'results', '-v7.3');
ewm_write_summary(results, cfg);

elapsed = toc(startTime);

% 写入可复现清单
manifestFile = fullfile(cfg.output.dir, 'reproducibility_manifest.txt');
fidManifest = fopen(manifestFile, 'w', 'n', 'UTF-8');
fprintf(fidManifest, '运行模式\t%s\n', mode);
fprintf(fidManifest, '运行时间\t%s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fidManifest, '模型路径\t%s\n', cfg.model.root);
fprintf(fidManifest, '输出路径\t%s\n', cfg.output.dir);
fprintf(fidManifest, '主模型网格\t%d x %d\n', model.nz, model.nx);
if isstruct(referenceComparison) && isfield(referenceComparison, 'referenceNz')
    fprintf(fidManifest, '参考解网格\t%d x %d\n', referenceComparison.referenceNz, referenceComparison.referenceNx);
else
    fprintf(fidManifest, '参考解网格\tN/A\n');
end
fprintf(fidManifest, 'dt\t%.12g\n', cfg.sim.dt);
fprintf(fidManifest, 'nt\t%d\n', cfg.sim.nt);
fprintf(fidManifest, '总模拟时长\t%.6f\n', cfg.sim.totalTime);
fprintf(fidManifest, '最大快照时刻\t%.6f\n', cfg.sim.maxSnapshotTime);
fprintf(fidManifest, 'PML层数\t%d\n', cfg.sim.nPml);
fprintf(fidManifest, 'Ricker主频\t%.2f\n', cfg.sim.f0);
fprintf(fidManifest, '震源位置\tz=%.2f m, xFraction=%.4f\n', ...
    cfg.sim.sourceDepthM, cfg.sim.sourceXFraction);
fprintf(fidManifest, 'Taylor系数\t%s\n', mat2str(standardCoeff, 8));
fprintf(fidManifest, '优化系数\t%s\n', mat2str(optimizedCoeff, 8));
fprintf(fidManifest, '优化最大频散误差\t%.6e\n', dispersion.optimized.maxAbsError);
fprintf(fidManifest, '目标误差阈值\t%.6e\n', cfg.coeff.targetError);
fprintf(fidManifest, '随机种子\t%d\n', cfg.coeff.sa.seed);
fprintf(fidManifest, 'MATLAB版本\t%s\n', version);
fprintf(fidManifest, '总耗时_秒\t%.2f\n', elapsed);
fclose(fidManifest);

% 输出文件存在性检查
expectedFiles = { ...
    'figures/marmousi_model.png', ...
    'figures/ricker_wavelet_time_spectrum.png', ...
    'figures/experiment_parameters_table.png', ...
    'experiment_parameters_table.csv', ...
    'figures/reference_solution_settings_table.png', ...
    'reference_solution_settings_table.csv', ...
    'figures/cfl_stability_table.png', ...
    'cfl_stability_table.csv', ...
    'figures/dispersion_curves.png', ...
    'figures/exp1_regular_vs_staggered.png', ...
    'figures/exp2_noabsorb_vs_pml.png', ...
    'figures/exp2_energy_noabsorb_vs_pml.png', ...
    'figures/exp2_pml_boundary_energy_time.png', ...
    'exp2_pml_boundary_energy_time.csv', ...
    'figures/exp3_standard_vs_minimax.png', ...
    'figures/exp3_taylor_optimized_difference.png', ...
    'exp3_taylor_optimized_difference_metrics.csv', ...
    'figures/exp3_taylor_optimized_difference_zoom.png', ...
    'exp3_taylor_optimized_difference_zoom.txt', ...
    'coefficients.txt', ...
    'summary.txt', ...
    'summary.mat', ...
    'reproducibility_manifest.txt' ...
};
if cfg.boundaryReference.enabled
    expectedFiles = [expectedFiles, { ...
        'figures/exp2_boundary_reference_error.png', ...
        'exp2_boundary_reference_metrics.txt' ...
    }];
end
if cfg.reference.enabled
    expectedFiles = [expectedFiles, { ...
        'figures/exp3_reference_error_comparison.png', ...
        'exp3_reference_error_metrics.txt', ...
        'exp3_snapshot_time_alignment.csv', ...
        'figures/sa_convergence.png' ...
    }];
end

missing = {};
for i = 1:numel(expectedFiles)
    if ~isfile(fullfile(cfg.output.dir, expectedFiles{i}))
        missing{end+1} = expectedFiles{i}; %#ok<AGROW>
    end
end

if ~isempty(missing)
    warning('ewm:MissingExpectedOutputs', ...
        '以下预期输出文件不存在：\n%s', strjoin(missing, newline));
end

fprintf('\n基于最大范数目标函数的优化差分系数（逐项）：\n');
for m = 1:numel(optimizedCoeff)
    fprintf('  c_%d = %.15g\n', m, optimizedCoeff(m));
end
fprintf('完整系数向量 = [');
fprintf(' %.15g', optimizedCoeff);
fprintf(' ]\n');

fprintf('\n完成。结果已写入：%s\n', cfg.output.dir);
fprintf('总耗时：%.2f 秒（约 %.2f 分钟）\n', elapsed, elapsed/60);
end
