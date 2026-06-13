function cfg = ewm_default_config(projectDir, mode)
%EWM_DEFAULT_CONFIG 为全部实验构建统一配置结构体。

cfg = struct();
cfg.mode = lower(char(mode));
cfg.projectDir = projectDir;

cfg.model = struct();
% 速度模型目录：项目内的 models/，下载者将自己的 Marmousi 模型文件放入即可。
cfg.model.root = fullfile(projectDir, 'models');
cfg.model.cropX = [];
cfg.model.cropZ = [];
cfg.model.enforcePhysicalVs = true;
cfg.model.vsOffset = 800;   % S 波速度统一加的偏移量（m/s），0 = 不加

cfg.coeff = struct();
cfg.coeff.order = 5;
cfg.coeff.khMax = 0.60 * pi;
cfg.coeff.samples = 500;
cfg.coeff.evalSamples = 1200;
% targetError 是绝对误差 max|k_num*Delta - k*Delta| 的阈值，
% 与 Zhang & Yao (2013) Eq. 23 一致；阈值 1e-4 与原文 Tables 1-2 完全相同口径。
cfg.coeff.targetError = 1.0e-4;
cfg.coeff.sa = struct();
cfg.coeff.sa.seed = 20260425;
cfg.coeff.sa.initialTemperature = 1.0;
cfg.coeff.sa.minTemperature = 1.0e-4;
cfg.coeff.sa.coolingRate = 0.992;
% 切换到绝对误差目标后，目标函数量级在 kh = 0.64π 处约为相对误差量级的 2 倍，
% 因此相比原相对误差版本需要更深的 SA + polish 搜索才能稳定压到 1e-4 以下。
cfg.coeff.sa.repeatPerTemperature = 160;
cfg.coeff.sa.restartCount = 40;
cfg.coeff.sa.temperaturePower = 0.65;
cfg.coeff.sa.acceptanceScale = 8.0;
cfg.coeff.sa.acceptanceFloor = 0.15;
cfg.coeff.sa.minPerturbationFraction = 0.002;
cfg.coeff.sa.perturbationScale = [0.030, 0.014, 0.005, 0.0015];
cfg.coeff.sa.restartScale = 0.7 * cfg.coeff.sa.perturbationScale;
cfg.coeff.sa.polishTrigger = 2.0;
cfg.coeff.sa.polishTrialsPerScale = 6000;
cfg.coeff.sa.polishMaxPasses = 8;
cfg.coeff.sa.polishStepRatios = [1.0, 0.5, 0.2, 0.05];
cfg.coeff.sa.validationSamples = 2000;

cfg.sim = struct();
cfg.sim.cfl = 0.42;
cfg.sim.nt = 260;
cfg.sim.dt = [];
cfg.sim.nPml = 12;
cfg.sim.f0 = 4.0;
cfg.sim.sourceDelayCycles = 1.5;
cfg.sim.sourceAmplitude = 1.0e8;
cfg.sim.sourceDepthM = 600;
cfg.sim.sourceXFraction = 0.50;
cfg.sim.snapshotFractions = [0.35, 0.65, 0.95];
cfg.snapshots = struct();
cfg.snapshots.exp1Times = [1.00, 1.10, 1.20];
cfg.snapshots.exp2Times = 1.65;
cfg.snapshots.exp3Times = [0.65, 0.95, 1.25];
cfg.reference = struct();
cfg.reference.enabled = false;
cfg.reference.spacing = [];
cfg.boundaryReference = struct();
cfg.boundaryReference.enabled = false;
cfg.boundaryReference.extraGrid = 25;
cfg.sim.energyStride = 1;
cfg.sim.outputDir = [];

switch cfg.mode
    case 'preview'
        cfg.model.spacing = 80;
        cfg.reference.enabled = true;
        cfg.reference.spacing = 40;
        cfg.boundaryReference.enabled = true;
        cfg.boundaryReference.extraGrid = 25;
        cfg.sim.nt = 360;
        cfg.sim.nPml = 12;
        cfg.sim.f0 = 4.0;
        cfg.output.dir = fullfile(projectDir, 'results_preview');

    case 'standard'
        cfg.model.spacing = 10;
        cfg.reference.enabled = false;
        cfg.reference.spacing = [];
        cfg.boundaryReference.enabled = true;
        cfg.boundaryReference.extraGrid = 100;
        cfg.sim.nt = 1600;
        cfg.sim.nPml = 30;
        cfg.sim.f0 = 12.0;
        cfg.output.dir = fullfile(projectDir, 'results_standard');

    case 'highfreq'
        cfg.model.spacing = 10;
        cfg.reference.enabled = false;
        cfg.reference.spacing = [];
        cfg.boundaryReference.enabled = true;
        cfg.boundaryReference.extraGrid = 80;
        cfg.sim.nt = 4000;
        cfg.sim.nPml = 20;
        cfg.sim.f0 = 25.0;
        cfg.output.dir = fullfile(projectDir, 'results_highfreq');

    otherwise
        error('未知运行模式 "%s"，请使用 "preview"、"standard" 或 "highfreq"。', cfg.mode);
end

cfg.sim.outputDir = cfg.output.dir;
end
