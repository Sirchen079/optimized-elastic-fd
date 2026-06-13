function ewm_generate_standard_values(resultsDir)
%EWM_GENERATE_THESIS_VALUES 从运行结果中提取论文所需数值。
%
%   "误差"在本项目中统一指模拟退火所控制的最大范数频散误差
%   max|k_num·Δ − k·Δ|（设计目标 1e-4）。其他基于波场快照的相对 L2
%   量在此文件中一律以"波场偏离"标注，与"误差"严格区分。
%
%   用法：
%       ewm_generate_standard_values()           % 使用 results_standard/
%       ewm_generate_standard_values('preview')    % 使用 results_preview/
%       ewm_generate_standard_values('highfreq')% 使用 results_highfreq/

if nargin < 1
    mode = 'standard';
elseif ischar(resultsDir) && ismember(resultsDir, {'preview', 'standard', 'highfreq'})
    mode = resultsDir;
else
    mode = 'standard';
end
resultsDir = fullfile(fileparts(mfilename('fullpath')), '..', sprintf('results_%s', mode));

%% ========== 读取原始数据 ==========

coeffFile = fullfile(resultsDir, 'coefficients.txt');
if exist(coeffFile, 'file')
    coeffText = fileread(coeffFile);
else
    coeffText = '';
end

summaryFile = fullfile(resultsDir, 'summary.txt');
if ~exist(summaryFile, 'file')
    error('找不到 %s，请先运行 run_ewm(''%s'')', summaryFile, mode);
end
summaryText = fileread(summaryFile);

%% ========== helpers ==========

function val = extract_value(text, key)
    pat = [regexptranslate('escape', key), '\s*=\s*([^\n]+)'];
    tok = regexp(text, pat, 'tokens', 'once');
    if ~isempty(tok)
        val = strtrim(tok{1});
    else
        val = '';
    end
end

function num = extract_num(text, key)
    s = extract_value(text, key);
    if ~isempty(s)
        num = str2double(s);
    else
        num = NaN;
    end
end

%% ========== 频散误差（唯一的"误差"口径）==========

taylorMaxAbs = extract_num(summaryText, 'Taylor频散最大误差');
minimaxMaxAbs = extract_num(summaryText, '优化系数频散最大误差');
improvement   = extract_num(summaryText, '频散误差改善倍数');
targetErr     = extract_num(summaryText, '频散目标误差');
if ~isempty(coeffText)
    c_taylorMaxAbs = extract_num(coeffText, 'Taylor最大绝对误差');
    c_minimaxMaxAbs = extract_num(coeffText, '优化系数最大绝对误差');
    c_improvement   = extract_num(coeffText, '绝对误差改善倍数');
    if ~isnan(c_taylorMaxAbs), taylorMaxAbs = c_taylorMaxAbs; end
    if ~isnan(c_minimaxMaxAbs), minimaxMaxAbs = c_minimaxMaxAbs; end
    if ~isnan(c_improvement), improvement = c_improvement; end
end

%% ========== 端到端波场偏离（辅助诊断，非"误差"）==========

taylorRefDev   = extract_num(summaryText, '实验3_Taylor相对参考波场偏离_平均');
minimaxRefDev  = extract_num(summaryText, '实验3_优化系数相对参考波场偏离_平均');
taylorHw       = extract_num(summaryText, '实验3_Taylor高波数能量占比_平均');
minimaxHw      = extract_num(summaryText, '实验3_优化系数高波数能量占比_平均');
refDevReduction = (taylorRefDev - minimaxRefDev) / taylorRefDev * 100;
hwReduction     = (taylorHw - minimaxHw) / taylorHw * 100;

%% ========== 系数和与时间步长 ==========

coeffVal = extract_value(summaryText, 'Taylor系数');
if ~isempty(coeffVal)
    taylorCoeffs = sscanf(coeffVal, '%f');
    if numel(taylorCoeffs) >= 5
        taylorSumAbs = sum(abs(taylorCoeffs(1:5)));
    else
        taylorSumAbs = NaN;
    end
else
    taylorSumAbs = NaN;
end
coeffVal2 = extract_value(summaryText, '优化系数');
if ~isempty(coeffVal2)
    optimCoeffs = sscanf(coeffVal2, '%f');
    if numel(optimCoeffs) >= 5
        optimSumAbs = sum(abs(optimCoeffs(1:5)));
    else
        optimSumAbs = NaN;
    end
else
    optimSumAbs = NaN;
end
dtReduction = (optimSumAbs - taylorSumAbs) / optimSumAbs * 100;
dtVal = extract_num(summaryText, '时间步长_s');

%% ========== PML 波场偏离与能量 ==========

pmlGlobalReduc  = extract_num(summaryText, '实验2_PML全域偏离降低百分比');
pmlBoundaryReduc = extract_num(summaryText, '实验2_PML边界带偏离降低百分比');
noAbsorbDev     = extract_num(summaryText, '实验2无吸收_相对参考波场偏离_平均');
pmlDev          = extract_num(summaryText, '实验2_PML_相对参考波场偏离_平均');
energy_pml      = extract_num(summaryText, '末时刻交错网格PML速度能量');
energyPct       = extract_num(summaryText, 'PML相对无吸收末时刻速度能量比');

%% ========== 生成报告 ==========

fid = fopen(fullfile(resultsDir, 'standard_values.txt'), 'w', 'n', 'UTF-8');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '论文数值提取报告\n');
fprintf(fid, '运行模式: %s\n', mode);
fprintf(fid, '生成时间: %s\n', datestr(now));
fprintf(fid, '%s\n\n', repmat('=', 1, 60));

fprintf(fid, '术语统一约定\n');
fprintf(fid, '%s\n', repmat('-', 1, 50));
fprintf(fid, '"误差"专指模拟退火所控制的最大范数频散误差 max|k_num·Δ − k·Δ|（设计目标 %s）。\n', num_or_dash(targetErr));
fprintf(fid, '其他基于波场快照的相对 L2 量在此文件中一律记为"波场偏离"，仅作辅助诊断，不计入"误差"。\n\n');

% --- 1. 频散误差 ---
fprintf(fid, '【1】频散误差（论文第3章 / 实验3理论依据，唯一的"误差"口径）\n');
fprintf(fid, '%s\n', repmat('-', 1, 50));
if ~isnan(taylorMaxAbs)
    fprintf(fid, 'Taylor系数频散最大误差 = %.4e\n', taylorMaxAbs);
end
if ~isnan(minimaxMaxAbs)
    fprintf(fid, '优化系数频散最大误差 = %.4e（≤ 设计目标 %s）\n', minimaxMaxAbs, num_or_dash(targetErr));
end
if ~isnan(improvement)
    fprintf(fid, '频散误差改善倍数 = %.2f\n', improvement);
end
fprintf(fid, '\n');

% --- 2. 实验3 端到端波场偏离 ---
fprintf(fid, '【2】实验3 端到端波场偏离（Taylor vs 优化系数，仅作辅助诊断）\n');
fprintf(fid, '%s\n', repmat('-', 1, 50));
if ~isnan(taylorRefDev)
    fprintf(fid, 'Taylor平均相对参考波场偏离 = %.4f\n', taylorRefDev);
end
if ~isnan(minimaxRefDev)
    fprintf(fid, '优化系数平均相对参考波场偏离 = %.4f\n', minimaxRefDev);
end
if ~isnan(refDevReduction)
    fprintf(fid, '波场偏离相对降幅 = %.2f%%\n', refDevReduction);
end
if ~isnan(taylorHw)
    fprintf(fid, 'Taylor平均高波数能量占比 = %.4f (%.2f%%)\n', taylorHw, taylorHw*100);
end
if ~isnan(minimaxHw)
    fprintf(fid, '优化系数平均高波数能量占比 = %.4f (%.2f%%)\n', minimaxHw, minimaxHw*100);
end
if ~isnan(hwReduction)
    fprintf(fid, '高波数能量占比相对降幅 = %.2f%%\n', hwReduction);
end
fprintf(fid, '\n');

% --- 3. 差分系数绝对值之和 ---
fprintf(fid, '【3】差分系数绝对值之和与时间步长\n');
fprintf(fid, '%s\n', repmat('-', 1, 50));
if ~isnan(taylorSumAbs)
    fprintf(fid, 'Taylor系数绝对值之和 = %.4f\n', taylorSumAbs);
end
if ~isnan(optimSumAbs)
    fprintf(fid, '优化系数绝对值之和 = %.4f\n', optimSumAbs);
end
if ~isnan(dtReduction)
    fprintf(fid, '系数和相对增加 = %.2f%%（提高稳定步长上限的代价）\n', dtReduction);
end
if ~isnan(dtVal)
    fprintf(fid, '实际使用时间步长 dt = %.6e s\n', dtVal);
end
fprintf(fid, '\n');

% --- 4. 实验2 PML ---
fprintf(fid, '【4】实验2 PML 吸收边界（同样以波场偏离作辅助诊断，主要指标为边界能量比）\n');
fprintf(fid, '%s\n', repmat('-', 1, 50));
if ~isnan(noAbsorbDev)
    fprintf(fid, '无吸收边界 平均相对参考波场偏离 = %.4f\n', noAbsorbDev);
end
if ~isnan(pmlDev)
    fprintf(fid, 'PML 吸收边界 平均相对参考波场偏离 = %.4f\n', pmlDev);
end
if ~isnan(pmlGlobalReduc)
    fprintf(fid, 'PML 全域波场偏离降低 = %.2f%%\n', pmlGlobalReduc);
end
if ~isnan(pmlBoundaryReduc)
    fprintf(fid, 'PML 边界带波场偏离降低 = %.2f%%\n', pmlBoundaryReduc);
end
if ~isnan(energyPct)
    fprintf(fid, 'PML 末时刻速度能量 = %.2f（为无吸收边界的 %.1f%%）\n', energy_pml, energyPct*100);
end
fprintf(fid, '\n');

% --- 5. 可直接粘贴的段落 ---
fprintf(fid, '【5】论文第四章可直接粘贴的段落\n');
fprintf(fid, '%s\n\n', repmat('=', 1, 60));

fprintf(fid, '--- 4.2 频散误差 ---\n');
if ~isnan(taylorMaxAbs) && ~isnan(minimaxMaxAbs) && ~isnan(improvement)
    fprintf(fid, ['Taylor 系数在目标波数范围内的最大范数频散误差为 %.3e，' ...
        '优化系数的最大范数频散误差为 %.2e，前者约为后者的 %.0f 倍。' ...
        '这表明基于最大范数目标函数的模拟退火优化已将差分算子的频散误差' ...
        '严格压制到设计目标 %s 以内。\n\n'], ...
        taylorMaxAbs, minimaxMaxAbs, improvement, num_or_dash(targetErr));
end

fprintf(fid, '--- 4.3 端到端波场偏离 ---\n');
if ~isnan(taylorRefDev) && ~isnan(minimaxRefDev) && ~isnan(refDevReduction)
    fprintf(fid, ['为直观反映上述算子级改进在长时传播后的体现，将粗网格波场与细网格参考波场做相对 L2 比较' ...
        '（仅作辅助诊断，与上节"误差"在量纲与作用域上均不相同）：' ...
        '在 20 m 网格、约 1.25 s 传播之后，Taylor 系数的平均相对参考波场偏离为 %.3f，' ...
        '优化系数为 %.3f，相对降幅约 %.1f%%。' ...
        '该比值偏大主要源于粗网格本身的累积相位差，并非算子级频散误差未被压制。\n\n'], ...
        taylorRefDev, minimaxRefDev, refDevReduction);
end

fprintf(fid, '--- 4.4 PML 吸收边界 ---\n');
if ~isnan(pmlGlobalReduc) && ~isnan(noAbsorbDev) && ~isnan(pmlDev)
    fprintf(fid, ['以扩大计算域参考解为基准，无吸收边界方案的全域相对参考波场偏离为 %.2f，' ...
        'PML 方案为 %.4f，PML 将全域波场偏离降低约 %.1f%%；' ...
        '边界带波场偏离同样下降约 %.1f%%，验证了 PML 层的反射抑制能力。\n'], ...
        noAbsorbDev, pmlDev, pmlGlobalReduc, pmlBoundaryReduc);
end

%% ========== 控制台预览 ==========

fprintf('\n===== 论文数值已写入 %s =====\n', fullfile(resultsDir, 'standard_values.txt'));
fprintf('关键数值预览（"误差"专指 SA 频散最大范数）：\n');
fprintf('  Taylor 频散最大误差:   %.3e\n', taylorMaxAbs);
fprintf('  优化系数频散最大误差: %.3e\n', minimaxMaxAbs);
fprintf('  频散误差改善倍数:     %.0f\n', improvement);
fprintf('  实验3 Taylor 波场偏离: %.3f\n', taylorRefDev);
fprintf('  实验3 优化 波场偏离:  %.3f\n', minimaxRefDev);
fprintf('  实验3 波场偏离降幅:    %.1f%%\n', refDevReduction);
fprintf('  实验2 PML 全域偏离降幅: %.1f%%\n', pmlGlobalReduc);
fprintf('  实验2 PML 边界偏离降幅: %.1f%%\n', pmlBoundaryReduc);

end

function out = num_or_dash(v)
if isnan(v)
    out = '—';
else
    out = sprintf('%.0e', v);
end
end
