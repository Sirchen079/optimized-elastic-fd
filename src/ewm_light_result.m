function out = ewm_light_result(in)
%EWM_LIGHT_RESULT 从摘要结构体中移除完整快照数据以节省空间。

out = in;
if isfield(out, 'snapshots')
    out = rmfield(out, 'snapshots');
end
end
