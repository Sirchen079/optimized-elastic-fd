%PLOT_SA_ERROR_THESIS_EVIDENCE Generate standard figures for SA absolute error evidence.
%
% This script does not rerun the wavefield simulation. It reads the saved
% coefficients and exp3 wavefield snapshots from results_standard, then creates
% figures where every quantitative "error" is the same SA absolute error:
%     |k_num * Delta - k * Delta|
% with target threshold 1e-4.

projectDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectDir, 'src'));

outputs = ewm_plot_sa_error_evidence('standard');

fprintf('\nGenerated SA absolute-error evidence figures:\n');
names = fieldnames(outputs.figures);
for i = 1:numel(names)
    value = outputs.figures.(names{i});
    if ischar(value) && ~isempty(value)
        fprintf('  %s\n', value);
    end
end
fprintf('Summary file:\n  %s\n', outputs.summaryFile);
