function idx = ewm_pick_snapshot_indices(sim, dt)
%EWM_PICK_SNAPSHOT_INDICES 将请求的快照时刻转换为时间步索引。

if isfield(sim, 'snapshotTimes') && ~isempty(sim.snapshotTimes)
    idx = round(sim.snapshotTimes(:).' / dt) + 1;
else
    idx = round(1 + (sim.nt - 1) * sim.snapshotFractions(:).');
end

idx = unique(max(1, min(sim.nt, idx)));
end
