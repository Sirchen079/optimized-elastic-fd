function redraw_standard_figures()
%REDRAW_THESIS_FIGURES 重绘论文所用的 8 张图（仅重画，不重新做数值实验）。
%
% 处理目标（输出文件名保持不变）：
%   1. exp1_regular_vs_staggered.png    常规/交错网格波场快照对比（裁剪+放大字体）
%   2. exp2_noabsorb_vs_pml.png         PML 吸收边界波场快照对比（裁剪+放大字体）
%   3. exp3_standard_vs_minimax.png     Taylor/优化系数波场快照对比（裁剪+放大字体）
%   4. marmousi_model.png               Marmousi 模型 Vp/Vs/密度（放大字体）
%   5. vs_model.png                     修改后的 Vs 模型（放大字体）
%   6. sa_error_theory_panel.png        优化系数频散分析四联图（放大字体）
%   7. sa_convergence.png               模拟退火收敛过程（放大字体）
%   8. ricker_wavelet_time_spectrum.png 雷克子波时域与频谱（放大字体）
%
% 数据全部取自 results_standard 下已保存的 .mat / coefficients.txt / summary.mat。

projectDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectDir, 'src'));
ewm_apply_chinese_style();

% 统一放大坐标轴、标题、图例、色标等字号，便于论文中清晰阅读。
set(groot, 'defaultAxesFontSize', 18);
set(groot, 'defaultAxesTitleFontSizeMultiplier', 1.12);
set(groot, 'defaultAxesLabelFontSizeMultiplier', 1.0);
set(groot, 'defaultColorbarFontSize', 16);
set(groot, 'defaultLegendFontSize', 16);
set(groot, 'defaultTextFontSize', 18);
restoreGroot = onCleanup(@reset_groot_fontsizes);

cfg = ewm_default_config(projectDir, 'standard');
ewm_ensure_dir(cfg.output.dir);
cfg.output.figuresDir = fullfile(cfg.output.dir, 'figures');
ewm_ensure_dir(cfg.output.figuresDir);

if cfg.model.spacing == 10
    ewm_build_marmousi_10m_cache(cfg.model.root, false);
end

fprintf('正在加载 Marmousi 模型...\n');
model = ewm_load_marmousi(cfg.model);

% --- 图 4 / 图 5：模型图 ---
fprintf('[1/6] 重画 marmousi_model.png 与 vs_model.png ...\n');
ewm_plot_model(model, cfg);

% --- 图 8：雷克子波（依赖实际运行的 dt/nt，取自 summary.mat）---
fprintf('[2/6] 重画 ricker_wavelet_time_spectrum.png ...\n');
summaryFile = fullfile(cfg.output.dir, 'summary.mat');
hasSummary = isfile(summaryFile);
if hasSummary
    S = load(summaryFile);
    if isfield(S, 'results') && isfield(S.results, 'config') && isfield(S.results.config, 'sim')
        cfg.sim.dt = S.results.config.sim.dt;
        cfg.sim.nt = S.results.config.sim.nt;
        if isfield(S.results.config.sim, 'totalTime')
            cfg.sim.totalTime = S.results.config.sim.totalTime;
        end
    end
end
if isempty(cfg.sim.dt)
    standardCoeff = ewm_fd_coefficients('standard', cfg.coeff.order);
    cfg.sim.dt = ewm_stable_dt(model, standardCoeff, cfg.sim.cfl);
end
if cfg.sim.nt <= 0
    cfg.sim.nt = ceil(1.65 / cfg.sim.dt) + 1;
end
ewm_plot_ricker_wavelet(cfg);

% --- 图 1：常规 vs 交错网格（裁剪左右空白）---
fprintf('[3/6] 重画 exp1_regular_vs_staggered.png ...\n');
exp1Regular = load_result(cfg.output.dir, 'exp1_regular_noabsorb.mat');
exp1Staggered = load_result(cfg.output.dir, 'exp1_staggered_noabsorb_standard.mat');
ewm_plot_comparison(model, exp1Regular, exp1Staggered, ...
    '常规网格', '交错网格', ...
    fullfile(cfg.output.figuresDir, 'exp1_regular_vs_staggered.png'), 2);

% --- 图 2：无吸收 vs PML（裁剪左右空白）---
fprintf('[4/6] 重画 exp2_noabsorb_vs_pml.png ...\n');
exp2NoAbsorb = load_result(cfg.output.dir, 'exp2_staggered_noabsorb_standard.mat');
exp2Pml = load_result(cfg.output.dir, 'exp2_staggered_pml_standard.mat');
ewm_plot_pml_snapshot_pair(model, exp2NoAbsorb, exp2Pml, ...
    fullfile(cfg.output.figuresDir, 'exp2_noabsorb_vs_pml.png'));

% --- 图 3：Taylor vs 优化系数（裁剪左右空白）---
fprintf('[5/6] 重画 exp3_standard_vs_minimax.png ...\n');
exp3Standard = load_result(cfg.output.dir, 'exp3_staggered_pml_standard.mat');
exp3Optimized = load_result(cfg.output.dir, 'exp3_staggered_pml_minimax.mat');
ewm_plot_comparison(model, exp3Standard, exp3Optimized, ...
    'Taylor 系数', '基于最大范数目标函数的优化系数', ...
    fullfile(cfg.output.figuresDir, 'exp3_standard_vs_minimax.png'));

% --- 图 6：优化系数频散分析四联图 ---
% 注：ewm_plot_sa_error_evidence 还会顺带重写 sa_error_kspace_map.png 与
% exp3_sa_wavefield_difference_clean.png 两张本次不需要改动的图；先备份后还原，
% 确保本次只改动 sa_error_theory_panel.png。
fprintf('[6/6] 重画 sa_error_theory_panel.png ...\n');
collateral = {fullfile(cfg.output.figuresDir, 'sa_error_kspace_map.png'), ...
              fullfile(cfg.output.figuresDir, 'exp3_sa_wavefield_difference_clean.png')};
backups = backup_files(collateral);
cleanupRestore = onCleanup(@() restore_files(backups));
ewm_plot_sa_error_evidence('standard');
% sa_convergence 不在 sa_error_evidence 内，单独重画（见下）。

% --- 图 7：模拟退火收敛过程 ---
fprintf('追加重画 sa_convergence.png ...\n');
optInfo = [];
if hasSummary && isfield(S, 'results') && isfield(S.results, 'optimization')
    optInfo = S.results.optimization;
end
if isempty(optInfo) || ~isfield(optInfo, 'trace')
    fprintf('summary.mat 中未找到优化 trace，重新运行模拟退火以获取收敛轨迹...\n');
    [~, optInfo] = ewm_optimize_minimax_coeffs( ...
        cfg.coeff.order, cfg.coeff.khMax, cfg.coeff.samples, cfg.coeff);
end
ewm_plot_sa_convergence(optInfo, cfg.coeff.targetError, ...
    fullfile(cfg.output.figuresDir, 'sa_convergence.png'));

fprintf('\n全部 8 张图已重新绘制完成。\n');
end

function result = load_result(outDir, fileName)
filePath = fullfile(outDir, fileName);
if ~exist(filePath, 'file')
    error('redraw:MissingSavedResult', '缺少已保存结果文件：%s', filePath);
end
loaded = load(filePath);
result = loaded.result;
end

function backups = backup_files(files)
backups = struct('orig', {}, 'temp', {});
for i = 1:numel(files)
    if isfile(files{i})
        tmp = [files{i}, '.redrawbak'];
        copyfile(files{i}, tmp);
        backups(end+1) = struct('orig', files{i}, 'temp', tmp); %#ok<AGROW>
    end
end
end

function restore_files(backups)
for i = 1:numel(backups)
    if isfile(backups(i).temp)
        copyfile(backups(i).temp, backups(i).orig);
        delete(backups(i).temp);
    end
end
end

function reset_groot_fontsizes()
set(groot, 'defaultAxesFontSize', 'remove');
set(groot, 'defaultAxesTitleFontSizeMultiplier', 'remove');
set(groot, 'defaultAxesLabelFontSizeMultiplier', 'remove');
set(groot, 'defaultColorbarFontSize', 'remove');
set(groot, 'defaultLegendFontSize', 'remove');
set(groot, 'defaultTextFontSize', 'remove');
end
