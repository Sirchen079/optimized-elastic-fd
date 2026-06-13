function ewm_save_figure(fig, outFile)
%EWM_SAVE_FIGURE 可靠地保存包含中文文字的图形。

[outDir, ~, ext] = fileparts(outFile);
ewm_ensure_dir(outDir);

% 确保图形已完成渲染以加载字体。
% 不可见图形有时会跳过字体加载。
drawnow_fig(fig);

% 显式对图形中的所有对象设置中文字体。
% 仅依赖 groot 默认值不可靠，tiledlayout、不可见图形和
% OpenGL 渲染器可能回退到非中文字体导致显示问号。
fontName = ewm_apply_chinese_style();
apply_font_to_all_objects(fig, fontName);

% 字体覆盖后再次刷新渲染，确保渲染器生效。
drawnow_fig(fig);

switch lower(ext)
    case '.png'
        % 优先使用 painters 矢量渲染器（字体保真度高）。
        % 若图形过于复杂（如大尺寸 imagesc），painters 可能失败，
        % 此时回退到图形当前的渲染器。
        oldRenderer = fig.Renderer;
        try
            fig.Renderer = 'painters';
            drawnow_fig(fig);
            print(fig, outFile, '-dpng', '-r900', '-painters');
        catch
            fig.Renderer = oldRenderer;
            drawnow_fig(fig);
            print(fig, outFile, '-dpng', '-r900');
        end
    case {'.jpg', '.jpeg'}
        print(fig, outFile, '-djpeg', '-r900');
    case '.pdf'
        print(fig, outFile, '-dpdf', '-r900');
    case '.eps'
        print(fig, outFile, '-depsc2', '-r900');
    otherwise
        try
            exportgraphics(fig, outFile, 'Resolution', 900);
        catch
            saveas(fig, outFile);
        end
end
end

function drawnow_fig(fig)
% 强制刷新指定图形的渲染。
% refresh(fig) 是更新单个图形的正确方式；
% drawnow 不接受图形句柄作为第一个参数。
try
    refresh(fig);
catch
    oldVisible = get(fig, 'Visible');
    if strcmpi(oldVisible, 'off')
        set(fig, 'Visible', 'on');
        drawnow;
        set(fig, 'Visible', 'off');
    else
        drawnow;
    end
end
end

function apply_font_to_all_objects(fig, fontName)
if isempty(fontName)
    return;
end
% findall 配合 -property 查找所有包含 FontName 属性的对象。
try
    objs = findall(fig, '-property', 'FontName');
    set(objs, 'FontName', fontName);
catch
    % 兼容不支持 -property 的旧版 MATLAB。
    types = {'axes', 'text', 'legend', 'colorbar', 'uitable'};
    for t = 1:numel(types)
        try
            set(findall(fig, 'Type', types{t}), 'FontName', fontName);
        catch
        end
    end
    try
        set(findall(fig, 'Type', 'textbox'), 'FontName', fontName);
    catch
    end
end
end
