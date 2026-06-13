function regen_standard_figs_hires()
%REGEN_THESIS_FIGS_HIRES 以高 DPI 重新生成论文中 8 张关键图。
%
% 涉及图：
%   1. marmousi_model.png
%   2. vs_model.png
%   3. ricker_wavelet_time_spectrum.png
%   4. exp1_regular_vs_staggered.png
%   5. exp2_noabsorb_vs_pml.png
%   6. exp3_standard_vs_minimax.png
%   7. sa_error_theory_panel.png
%   8. sa_convergence.png
%
% 依赖 results_standard 下已保存的 .mat 结果和 coefficients.txt、summary.mat。
% 不重新跑数值实验，只用已有结果重画图。

projectDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectDir, 'src'));
ewm_apply_chinese_style();

cfg = ewm_default_config(projectDir, 'standard');
ewm_ensure_dir(cfg.output.dir);
cfg.output.figuresDir = fullfile(cfg.output.dir, 'figures');
ewm_ensure_dir(cfg.output.figuresDir);

if (cfg.reference.enabled && cfg.reference.spacing == 10) || cfg.model.spacing == 10
    ewm_build_marmousi_10m_cache(cfg.model.root, false);
end

fprintf('正在加载 Marmousi 模型...\n');
model = ewm_load_marmousi(cfg.model);

fprintf('[1/2] 正在重画 marmousi_model.png 与 vs_model.png ...\n');
ewm_plot_model(model, cfg);

fprintf('[3] 正在重画 ricker_wavelet_time_spectrum.png ...\n');
% ewm_plot_ricker_wavelet 依赖 cfg.sim.dt 与 cfg.sim.nt。
% 这里用 summary.mat 中的实际运行参数，确保和正式运行一致。
summaryFile = fullfile(cfg.output.dir, 'summary.mat');
hasSummary = isfile(summaryFile);
if hasSummary
    S = load(summaryFile);
    if isfield(S, 'results') && isfield(S.results, 'config') ...
            && isfield(S.results.config, 'sim')
        cfg.sim.dt = S.results.config.sim.dt;
        cfg.sim.nt = S.results.config.sim.nt;
        if isfield(S.results.config.sim, 'totalTime')
            cfg.sim.totalTime = S.results.config.sim.totalTime;
        end
    end
end
if isempty(cfg.sim.dt)
    % 兜底：用 Taylor 系数估算 dt。
    standardCoeff = ewm_fd_coefficients('standard', cfg.coeff.order);
    cfg.sim.dt = ewm_stable_dt(model, standardCoeff, cfg.sim.cfl);
end
if cfg.sim.nt <= 0
    cfg.sim.nt = ceil(1.65 / cfg.sim.dt) + 1;
end
ewm_plot_ricker_wavelet(cfg);

fprintf('[4] 正在重画 exp1_regular_vs_staggered.png ...\n');
regularNoAbsorb = load_result(cfg.output.dir, 'exp1_regular_noabsorb.mat');
staggeredNoAbsorbExp1 = load_result(cfg.output.dir, 'exp1_staggered_noabsorb_standard.mat');
ewm_plot_comparison(model, regularNoAbsorb, staggeredNoAbsorbExp1, ...
    '常规网格', '交错网格', ...
    fullfile(cfg.output.figuresDir, 'exp1_regular_vs_staggered.png'), 2);

fprintf('[5] 正在重画 exp2_noabsorb_vs_pml.png ...\n');
exp2NoAbsorb = load_result(cfg.output.dir, 'exp2_staggered_noabsorb_standard.mat');
exp2Pml = load_result(cfg.output.dir, 'exp2_staggered_pml_standard.mat');
ewm_plot_pml_snapshot_pair(model, exp2NoAbsorb, exp2Pml, ...
    fullfile(cfg.output.figuresDir, 'exp2_noabsorb_vs_pml.png'));

fprintf('[6] 正在重画 exp3_standard_vs_minimax.png ...\n');
exp3Standard = load_result(cfg.output.dir, 'exp3_staggered_pml_standard.mat');
exp3Optimized = load_result(cfg.output.dir, 'exp3_staggered_pml_minimax.mat');
ewm_plot_comparison(model, exp3Standard, exp3Optimized, ...
    'Taylor 系数', '基于最大范数目标函数的优化系数', ...
    fullfile(cfg.output.figuresDir, 'exp3_standard_vs_minimax.png'));

fprintf('[7] 正在重画 sa_error_theory_panel.png 及其同组图 ...\n');
ewm_plot_sa_error_evidence('standard');

fprintf('[8] 正在重画 sa_convergence.png ...\n');
optInfo = [];
if hasSummary && isfield(S, 'results') && isfield(S.results, 'optimization')
    optInfo = S.results.optimization;
end
if isempty(optInfo) || ~isfield(optInfo, 'trace')
    fprintf('summary.mat 中未找到优化 trace，正在重新运行模拟退火以获取收敛轨迹（可能需要数秒）...\n');
    [~, optInfo] = ewm_optimize_minimax_coeffs( ...
        cfg.coeff.order, cfg.coeff.khMax, cfg.coeff.samples, cfg.coeff);
end
ewm_plot_sa_convergence(optInfo, cfg.coeff.targetError, ...
    fullfile(cfg.output.figuresDir, 'sa_convergence.png'));

fprintf('\n所有 8 张图已以高 DPI 重新生成完成。\n');
end

function result = load_result(outDir, fileName)
filePath = fullfile(outDir, fileName);
if ~exist(filePath, 'file')
    error('regen:MissingSavedResult', '缺少已保存结果文件：%s', filePath);
end
loaded = load(filePath);
result = loaded.result;
end
