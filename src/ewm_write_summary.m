function ewm_write_summary(results, cfg)
%EWM_WRITE_SUMMARY 写入论文用文本摘要报告。

outFile = fullfile(cfg.output.dir, 'summary.txt');
fid = fopen(outFile, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'ewm 弹性波数值模拟结果摘要\n');
fprintf(fid, '运行模式 = %s\n', mode_text(cfg.mode));
fprintf(fid, '模型网格点数_z = %d\n', results.model.nz);
fprintf(fid, '模型网格点数_x = %d\n', results.model.nx);
fprintf(fid, 'dz = %.6g\n', results.model.dz);
fprintf(fid, 'dx = %.6g\n', results.model.dx);
fprintf(fid, '时间步长_s = %.15g\n', cfg.sim.dt);
fprintf(fid, '时间步数 = %d\n', cfg.sim.nt);
fprintf(fid, 'Ricker主频_Hz = %.6g\n\n', cfg.sim.f0);
fprintf(fid, '实验1快照时刻_s =');
fprintf(fid, ' %.6g', cfg.snapshots.exp1Times);
fprintf(fid, '\n实验2快照时刻_s =');
fprintf(fid, ' %.6g', cfg.snapshots.exp2Times);
fprintf(fid, '\n实验3快照时刻_s =');
fprintf(fid, ' %.6g', cfg.snapshots.exp3Times);
fprintf(fid, '\n\n');

fprintf(fid, 'Taylor系数 =');
fprintf(fid, ' %.15g', results.coefficients.standard);
fprintf(fid, '\n优化系数 =');
fprintf(fid, ' %.15g', results.coefficients.optimized);
fprintf(fid, '\n\n');
fprintf(fid, '优化系数逐项\n');
for m = 1:numel(results.coefficients.optimized)
    fprintf(fid, 'c_%d = %.15g\n', m, results.coefficients.optimized(m));
end
fprintf(fid, '\n');

if isfield(results, 'optimization') && isfield(results.optimization, 'method')
    fprintf(fid, '优化方法 = %s\n', optimizer_text(results.optimization.method));
    fprintf(fid, '优化随机种子 = %.15g\n', results.optimization.seed);
    fprintf(fid, '优化验证最大误差 = %.15g\n\n', ...
        results.optimization.validationMaxError);
end

fprintf(fid, '频散最大kh除以pi = %.15g\n', cfg.coeff.khMax / pi);
fprintf(fid, '频散目标误差 = %.15g\n', cfg.coeff.targetError);
fprintf(fid, '常规网格频散最大误差 = %.15g\n', ...
    results.dispersion.regular.maxAbsError);
fprintf(fid, 'Taylor频散最大误差 = %.15g\n', ...
    results.dispersion.standard.maxAbsError);
fprintf(fid, '优化系数频散最大误差 = %.15g\n', ...
    results.dispersion.optimized.maxAbsError);
fprintf(fid, '频散误差改善倍数 = %.15g\n\n', ...
    results.dispersion.improvement);

fprintf(fid, '末时刻常规网格无吸收边界能量占比 = %.15g\n', ...
    results.regularNoAbsorb.energy.boundaryRatio(end));
fprintf(fid, '末时刻交错网格无吸收边界能量占比 = %.15g\n', ...
    results.staggeredNoAbsorb.energy.boundaryRatio(end));
fprintf(fid, '末时刻交错网格PML边界能量占比 = %.15g\n', ...
    results.staggeredPmlStandard.energy.boundaryRatio(end));
fprintf(fid, '末时刻优化系数PML边界能量占比 = %.15g\n', ...
    results.staggeredPmlOptimized.energy.boundaryRatio(end));
fprintf(fid, '\n');
fprintf(fid, '末时刻常规网格无吸收速度能量 = %.15g\n', ...
    results.regularNoAbsorb.energy.totalVelocity(end));
fprintf(fid, '末时刻交错网格无吸收速度能量 = %.15g\n', ...
    results.staggeredNoAbsorb.energy.totalVelocity(end));
fprintf(fid, '末时刻交错网格PML速度能量 = %.15g\n', ...
    results.staggeredPmlStandard.energy.totalVelocity(end));
fprintf(fid, '末时刻优化系数PML速度能量 = %.15g\n', ...
    results.staggeredPmlOptimized.energy.totalVelocity(end));
fprintf(fid, 'PML相对无吸收末时刻速度能量比 = %.15g\n', ...
    results.staggeredPmlStandard.energy.totalVelocity(end) / ...
    (results.staggeredNoAbsorb.energy.totalVelocity(end) + eps));

fprintf(fid, '\n论文写作注意事项\n');
fprintf(fid, '----------------\n');
fprintf(fid, '1. 时间步长说明\n');
fprintf(fid, '   主模型模拟采用统一时间步长 dt = %.15g s，该值取 Taylor 系数稳定时间步长和基于最大范数目标函数的优化系数稳定时间步长中的较小者，因此两类系数对比实验具有相同时间离散条件。\n', cfg.sim.dt);
fprintf(fid, '\n');
fprintf(fid, '2. PML外扩说明\n');
fprintf(fid, '   PML吸收层通过在有效物理模型外侧扩展网格实现，波场输出时裁剪回原始物理模型区域。因此震源位于有效模型内部，而不是位于PML区域内。\n');
fprintf(fid, '\n');
fprintf(fid, '3. Vs = 0 说明\n');
fprintf(fid, '   Marmousi模型局部区域存在 Vs = 0 的近似流体或极低剪切模量区域，参数转换时采用 μ = ρVs²。当 Vs = 0 时，对应区域剪切模量 μ = 0。稳定性主要由最大 P 波速度控制。\n');
fprintf(fid, '\n');
fprintf(fid, '4. 误差定义统一约定\n');
fprintf(fid, '   全文"误差"专指模拟退火所控制的最大范数频散误差 max|k_num·Δ − k·Δ|（设计目标 1e-4），即 Taylor / 优化系数对一阶空间差分算子的频域偏差。\n');
fprintf(fid, '\n');
fprintf(fid, '5. 常规网格对比说明\n');
fprintf(fid, '   常规网格与交错网格对比实验仅作为基线现象展示，不作为严格同阶精度评估。后续所有核心实验均采用交错网格有限差分框架。\n');
end

function text = optimizer_text(value)
if isstring(value)
    value = char(value);
end
if strcmp(value, 'simulated_annealing_maximum_norm')
    text = '最大范数目标模拟退火';
else
    text = value;
end
end

function text = mode_text(value)
if isstring(value)
    value = char(value);
end
if strcmp(value, 'standard')
    text = '论文模式';
elseif strcmp(value, 'preview')
    text = '快速验证模式';
else
    text = value;
end
end
