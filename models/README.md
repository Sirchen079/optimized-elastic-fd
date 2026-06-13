# 速度模型目录

本目录存放弹性波模拟所需的速度模型文件。代码通过相对路径定位到 `models/`，**下载项目后将自己的模型文件放入此处即可运行，无需修改任何代码路径**。

## 需要放入的文件

本项目默认使用 **AGL 弹性 Marmousi 模型**。按运行模式不同，需要的文件如下：

| 运行模式 | 网格间距 | 需要的文件 |
|----------|----------|------------|
| `preview` | 80 m | `marmousi_cache_20m.mat`（代码自动抽稀到 80 m） |
| `standard` | 10 m | `marmousi_cache_10m.mat` |

### 方式一：直接放缓存文件（推荐）

将预先生成的 `.mat` 缓存文件放入本目录：
- `marmousi_cache_20m.mat` — 基础缓存（`preview` 模式必需）
- `marmousi_cache_10m.mat` — 细网格缓存（`standard` 模式必需）

缓存文件为 MATLAB `-v7.3` 格式，需包含字段：`vp`、`vs`、`rho`（P 波速度、S 波速度、密度，单位 m/s）、`dx`、`dz`（网格间距，单位 m）。

### 方式二：从原始 SEG-Y 自动构建

若放入原始 1.25 m SEG-Y 文件，首次以 `standard` 模式运行时代码会自动重采样生成 10 m 缓存：
- `MODEL_P-WAVE_VELOCITY_1.25m.segy`
- `MODEL_S-WAVE_VELOCITY_1.25m.segy`
- `MODEL_DENSITY_1.25m.segy`

## 获取 Marmousi 模型

AGL 弹性 Marmousi 模型可从以下来源获取：
- Agile Geoscience Marmousi2：https://github.com/agile-geoscience/marmousi2

也可使用任意标准 Marmousi 2 速度模型，按下方格式整理为 `.mat` 缓存。

## 使用自定义模型

若使用自己的速度模型，按缓存格式生成 `.mat` 放入本目录即可：

```matlab
vp  = ...;  % P 波速度，nz × nx 矩阵，单位 m/s
vs  = ...;  % S 波速度，nz × nx 矩阵，单位 m/s
rho = ...;  % 密度，    nz × nx 矩阵，单位 kg/m³
dx  = 20;   % 水平网格间距 (m)
dz  = 20;   % 垂直网格间距 (m)
save('marmousi_cache_20m.mat', 'vp','vs','rho','dx','dz', '-v7.3');
```

> 文件名中的间距数值必须与实际 `dx` / `dz` 一致，代码据此选择对应缓存。
