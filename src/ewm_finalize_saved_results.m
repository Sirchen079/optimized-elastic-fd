function results = ewm_finalize_saved_results(mode)
%EWM_FINALIZE_SAVED_RESULTS 从已保存的 MAT 文件重新生成图表和摘要。

if nargin < 1 || isempty(mode)
    mode = 'standard';
end

projectDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectDir, 'src'));
ewm_apply_chinese_style();

cfg = ewm_default_config(projectDir, mode);
if cfg.reference.enabled && cfg.reference.spacing == 10
    ewm_build_marmousi_10m_cache(cfg.model.root, false);
end
ewm_ensure_dir(cfg.output.dir);
cfg.output.figuresDir = fullfile(cfg.output.dir, 'figures');
ewm_ensure_dir(cfg.output.figuresDir);
startTime = tic;

model = ewm_load_marmousi(cfg.model);
standardCoeff = ewm_fd_coefficients('standard', cfg.coeff.order);
[optimizedCoeff, optInfo] = ewm_optimize_minimax_coeffs( ...
    cfg.coeff.order, cfg.coeff.khMax, cfg.coeff.samples, cfg.coeff);
cfg.sim.dtTaylorLimit = ewm_stable_dt(model, standardCoeff, cfg.sim.cfl);
cfg.sim.dtOptimizedLimit = ewm_stable_dt(model, optimizedCoeff, cfg.sim.cfl);
cfg.sim.dt = min(cfg.sim.dtTaylorLimit, cfg.sim.dtOptimizedLimit);
allSnapshotTimes = [ ...
    cfg.snapshots.exp1Times(:); ...
    cfg.snapshots.exp2Times(:); ...
    cfg.snapshots.exp3Times(:)];
cfg.sim.maxSnapshotTime = max(allSnapshotTimes);
cfg.sim.nt = max(cfg.sim.nt, ceil(cfg.sim.maxSnapshotTime / cfg.sim.dt) + 1);
cfg.sim.totalTime = cfg.sim.dt * (cfg.sim.nt - 1);

dispersion = ewm_dispersion_metrics(standardCoeff, optimizedCoeff, ...
    cfg.coeff.khMax, cfg.coeff.evalSamples);
ewm_plot_model(model, cfg);
ewm_plot_ricker_wavelet(cfg);
ewm_plot_dispersion(standardCoeff, optimizedCoeff, dispersion, cfg, optInfo);

regularNoAbsorb = load_result(cfg.output.dir, 'exp1_regular_noabsorb.mat');
staggeredNoAbsorbExp1 = load_result(cfg.output.dir, 'exp1_staggered_noabsorb_standard.mat');
staggeredNoAbsorbExp2 = load_result(cfg.output.dir, 'exp2_staggered_noabsorb_standard.mat');
staggeredPmlStandardExp2 = load_result(cfg.output.dir, 'exp2_staggered_pml_standard.mat');
boundaryReferenceResult = load_result(cfg.output.dir, 'exp2_boundary_large_domain_reference.mat');
staggeredPmlStandardExp3 = load_result(cfg.output.dir, 'exp3_staggered_pml_standard.mat');
staggeredPmlOptimized = load_result(cfg.output.dir, 'exp3_staggered_pml_minimax.mat');
referenceResult = load_result(cfg.output.dir, 'exp3_reference_finer_grid.mat');

simBase = cfg.sim;
simBase.standardCoeff = standardCoeff;
simBase.optimizedCoeff = optimizedCoeff;

ewm_plot_comparison(model, regularNoAbsorb, staggeredNoAbsorbExp1, ...
    '常规网格', '交错网格', ...
    fullfile(cfg.output.figuresDir, 'exp1_regular_vs_staggered.png'), 2);
ewm_plot_pml_snapshot_pair(model, staggeredNoAbsorbExp2, staggeredPmlStandardExp2, ...
    fullfile(cfg.output.figuresDir, 'exp2_noabsorb_vs_pml.png'));
ewm_plot_energy(staggeredNoAbsorbExp2, staggeredPmlStandardExp2, ...
    '无吸收边界', 'PML 吸收边界', ...
    fullfile(cfg.output.figuresDir, 'exp2_energy_noabsorb_vs_pml.png'));
ewm_plot_pml_boundary_energy(staggeredNoAbsorbExp2, staggeredPmlStandardExp2, ...
    fullfile(cfg.output.figuresDir, 'exp2_pml_boundary_energy_time.png'), ...
    fullfile(cfg.output.dir, 'exp2_pml_boundary_energy_time.csv'));

simExp2 = simBase;
simExp2.snapshotTimes = cfg.snapshots.exp2Times;
simExp2.nt = max(simExp2.nt, ceil(max(simExp2.snapshotTimes) / simExp2.dt) + 1);
[modelBoundaryRef, cropInfo] = ewm_make_boundary_reference_case( ...
    model, simExp2, cfg.boundaryReference.extraGrid);
boundaryReferenceComparison = ewm_plot_boundary_reference_error(model, ...
    staggeredNoAbsorbExp2, staggeredPmlStandardExp2, modelBoundaryRef, ...
    boundaryReferenceResult, cropInfo, ...
    fullfile(cfg.output.figuresDir, 'exp2_boundary_reference_error.png'), ...
    fullfile(cfg.output.dir, 'exp2_boundary_reference_metrics.txt'));

ewm_plot_comparison(model, staggeredPmlStandardExp3, staggeredPmlOptimized, ...
    'Taylor 系数', '基于最大范数目标函数的优化系数', fullfile(cfg.output.figuresDir, 'exp3_standard_vs_minimax.png'));
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

modelRef = [];
referenceComparison = [];
if cfg.reference.enabled && cfg.reference.spacing < cfg.model.spacing
    refModelCfg = cfg.model;
    refModelCfg.spacing = cfg.reference.spacing;
    modelRef = ewm_load_marmousi(refModelCfg);
    simExp3 = simBase;
    simExp3.snapshotTimes = cfg.snapshots.exp3Times;
    simRef = simExp3;
    simRef.dt = ewm_stable_dt(modelRef, standardCoeff, cfg.sim.cfl);
    totalTime = cfg.sim.dt * (cfg.sim.nt - 1);
    simRef.nt = ceil(totalTime / simRef.dt) + 1;
    simRef.nPml = max(simExp3.nPml, round(simExp3.nPml * model.dx / modelRef.dx));
    referenceComparison = ewm_plot_reference_error(model, ...
        staggeredPmlStandardExp3, staggeredPmlOptimized, modelRef, referenceResult, ...
        fullfile(cfg.output.figuresDir, 'exp3_reference_error_comparison.png'), ...
        fullfile(cfg.output.dir, 'exp3_reference_error_metrics.txt'));
    referenceComparison = enrich_reference_metrics(referenceComparison, modelRef, simRef);
end

requiredTables = ewm_write_required_tables(model, cfg, standardCoeff, ...
    optimizedCoeff, optInfo, dispersion, boundaryReferenceComparison, referenceComparison);

results = build_results_struct(cfg, model, standardCoeff, optimizedCoeff, optInfo, ...
    dispersion, regularNoAbsorb, staggeredNoAbsorbExp2, staggeredPmlStandardExp2, ...
    boundaryReferenceComparison, staggeredPmlStandardExp3, staggeredPmlOptimized, ...
    referenceComparison, exp3WavefieldDifference, exp3ZoomInfo, requiredTables);

summaryFile = fullfile(cfg.output.dir, 'summary.mat');
save(summaryFile, 'results', '-v7.3');
elapsed = toc(startTime);
fprintf('总耗时：%.2f 秒（约 %.2f 分钟）\n', elapsed, elapsed/60);
ewm_write_summary(results, cfg);
end

function result = load_result(outDir, fileName)
filePath = fullfile(outDir, fileName);
if ~exist(filePath, 'file')
    error('ewm:MissingSavedResult', '缺少已保存结果文件：%s', filePath);
end
loaded = load(filePath);
result = loaded.result;
end

function metrics = enrich_reference_metrics(metrics, modelRef, simRef)
metrics.referenceDt = simRef.dt;
metrics.referenceNt = simRef.nt;
metrics.referencePmlLayers = simRef.nPml;
metrics.referenceSource = modelRef.source;
metrics.referenceCacheSpacing = modelRef.cacheSpacing;
metrics.referenceMaxVp = max(modelRef.vp(:));
if isfield(modelRef, 'sourceSpacing')
    metrics.referenceSourceSpacing = modelRef.sourceSpacing;
end
if isfield(modelRef, 'resample')
    metrics.referenceResample = modelRef.resample;
end
end

function results = build_results_struct(cfg, model, standardCoeff, optimizedCoeff, optInfo, ...
    dispersion, regularNoAbsorb, staggeredNoAbsorbExp2, staggeredPmlStandardExp2, ...
    boundaryReferenceComparison, staggeredPmlStandardExp3, staggeredPmlOptimized, ...
    referenceComparison, exp3WavefieldDifference, exp3ZoomInfo, requiredTables)
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
end
