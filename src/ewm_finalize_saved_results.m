function results = ewm_finalize_saved_results(mode)
%EWM_FINALIZE_SAVED_RESULTS 从已保存的 MAT 文件重新生成图表和摘要。

if nargin < 1 || isempty(mode)
    mode = 'standard';
end

projectDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectDir, 'src'));
ewm_apply_chinese_style();

cfg = ewm_default_config(projectDir, mode);
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
staggeredPmlStandardExp3 = load_result(cfg.output.dir, 'exp3_staggered_pml_standard.mat');
staggeredPmlOptimized = load_result(cfg.output.dir, 'exp3_staggered_pml_minimax.mat');

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

requiredTables = ewm_write_required_tables(model, cfg, standardCoeff, ...
    optimizedCoeff, optInfo, dispersion);

results = build_results_struct(cfg, model, standardCoeff, optimizedCoeff, optInfo, ...
    dispersion, regularNoAbsorb, staggeredNoAbsorbExp2, staggeredPmlStandardExp2, ...
    staggeredPmlStandardExp3, staggeredPmlOptimized, ...
    exp3WavefieldDifference, exp3ZoomInfo, requiredTables);

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

function results = build_results_struct(cfg, model, standardCoeff, optimizedCoeff, optInfo, ...
    dispersion, regularNoAbsorb, staggeredNoAbsorbExp2, staggeredPmlStandardExp2, ...
    staggeredPmlStandardExp3, staggeredPmlOptimized, ...
    exp3WavefieldDifference, exp3ZoomInfo, requiredTables)
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
results.staggeredPmlStandardCoeffExp3 = ewm_light_result(staggeredPmlStandardExp3);
results.staggeredPmlOptimized = ewm_light_result(staggeredPmlOptimized);
results.exp3WavefieldDifference = exp3WavefieldDifference;
results.exp3ZoomInfo = exp3ZoomInfo;
results.requiredTables = requiredTables;
end
