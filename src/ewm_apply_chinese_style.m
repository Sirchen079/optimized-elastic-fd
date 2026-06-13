function fontName = ewm_apply_chinese_style()
%EWM_APPLY_CHINESE_STYLE 检测并设置支持中文的字体。
% 返回所选字体名称（未找到则为空）。

persistent cachedFontName
if ~isempty(cachedFontName)
    fontName = cachedFontName;
    apply_to_groot(fontName);
    return;
end

availableFonts = listfonts;

% 按 Windows/MATLAB 常见安装情况排列优先级。
% 同时包含英文名和中文名，因为 listfonts 返回的名称
% 取决于系统区域设置和 MATLAB 版本。
candidates = {'Microsoft YaHei', 'SimHei', 'SimSun', 'NSimSun', ...
              'FangSong', 'KaiTi', 'STHeiti', 'STSong', ...
              'Noto Sans CJK SC', 'Source Han Sans SC', ...
              'WenQuanYi Micro Hei', 'WenQuanYi Zen Hei', ...
              'Heiti SC', 'PingFang SC', 'Arial Unicode MS', ...
              '微软雅黑', '黑体', '宋体', '新宋体', ...
              '仿宋', '楷体', '华文黑体', '华文宋体'};

fontName = '';
for i = 1:numel(candidates)
    idx = find(strcmpi(candidates{i}, availableFonts), 1);
    if ~isempty(idx)
        % 使用 listfonts 返回的原始拼写，避免大小写或空格不匹配。
        fontName = availableFonts{idx};
        break;
    end
end

if isempty(fontName)
    warning('ewm:apply_chinese_style', ...
            'No Chinese font found among candidates. Chinese text may display as boxes or question marks.');
    cachedFontName = '';
    return;
end

cachedFontName = fontName;
apply_to_groot(fontName);
end

function apply_to_groot(fontName)
if isempty(fontName)
    return;
end
set(groot, 'defaultAxesFontName', fontName);
set(groot, 'defaultTextFontName', fontName);
set(groot, 'defaultLegendFontName', fontName);
set(groot, 'defaultColorbarFontName', fontName);
end
