function regen_two_figs()
%REGEN_TWO_FIGS 仅重新生成 marmousi_model.png 和 exp1_regular_vs_staggered.png。

projectDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectDir, 'src'));
ewm_apply_chinese_style();

cfg = ewm_default_config(projectDir, 'standard');
if cfg.reference.enabled && cfg.reference.spacing == 10
    ewm_build_marmousi_10m_cache(cfg.model.root, false);
end
ewm_ensure_dir(cfg.output.dir);
cfg.output.figuresDir = fullfile(cfg.output.dir, 'figures');
ewm_ensure_dir(cfg.output.figuresDir);

model = ewm_load_marmousi(cfg.model);
ewm_plot_model(model, cfg);
fprintf('Regenerated: %s\n', fullfile(cfg.output.figuresDir, 'marmousi_model.png'));

regularNoAbsorb = load_result(cfg.output.dir, 'exp1_regular_noabsorb.mat');
staggeredNoAbsorbExp1 = load_result(cfg.output.dir, 'exp1_staggered_noabsorb_standard.mat');
ewm_plot_comparison(model, regularNoAbsorb, staggeredNoAbsorbExp1, ...
    '常规网格', '交错网格', ...
    fullfile(cfg.output.figuresDir, 'exp1_regular_vs_staggered.png'), 2);
fprintf('Regenerated: %s\n', fullfile(cfg.output.figuresDir, 'exp1_regular_vs_staggered.png'));

staggeredPmlStandardExp3 = load_result(cfg.output.dir, 'exp3_staggered_pml_standard.mat');
staggeredPmlOptimized = load_result(cfg.output.dir, 'exp3_staggered_pml_minimax.mat');
ewm_plot_comparison(model, staggeredPmlStandardExp3, staggeredPmlOptimized, ...
    'Taylor 系数', '基于最大范数目标函数的优化系数', ...
    fullfile(cfg.output.figuresDir, 'exp3_standard_vs_minimax.png'));
fprintf('Regenerated: %s\n', fullfile(cfg.output.figuresDir, 'exp3_standard_vs_minimax.png'));
end

function result = load_result(outDir, fileName)
filePath = fullfile(outDir, fileName);
if ~exist(filePath, 'file')
    error('regen:MissingSavedResult', '缺少已保存结果文件：%s', filePath);
end
loaded = load(filePath);
result = loaded.result;
end
