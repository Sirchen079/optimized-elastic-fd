function ewm_plot_ricker_wavelet(cfg)
%EWM_PLOT_RICKER_WAVELET 保存 Ricker 子波时域波形与频谱图。

ewm_apply_chinese_style();

t = (0:cfg.sim.nt-1) * cfg.sim.dt;
wavelet = ewm_ricker(t, cfg.sim.f0, cfg.sim.sourceDelayCycles);
source = cfg.sim.sourceAmplitude * wavelet;

nfft = 2 ^ nextpow2(max(numel(wavelet), 1024));
freq = (0:nfft/2).' / (nfft * cfg.sim.dt);
spectrum = abs(fft(wavelet, nfft));
spectrum = spectrum(1:nfft/2+1).';
spectrum = spectrum(:);
if max(spectrum) > 0
    spectrum = spectrum / max(spectrum);
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1180, 460]);
tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(t, wavelet, 'k-', 'LineWidth', 1.8);
grid on;
box on;
xlabel('时间 (s)');
ylabel('归一化振幅');
title(sprintf('Ricker 子波，主频 = %.3g Hz', cfg.sim.f0));
set(gca, 'LineWidth', 0.9, 'GridAlpha', 0.18);

nexttile;
plot(freq, spectrum, 'Color', [0.05, 0.32, 0.62], 'LineWidth', 1.8);
grid on;
box on;
xlabel('频率 (Hz)');
ylabel('归一化振幅谱');
title('Ricker 子波振幅谱');
set(gca, 'LineWidth', 0.9, 'GridAlpha', 0.18);
xlim([0, max(cfg.sim.f0 * 4, cfg.sim.f0 + eps)]);

ewm_save_figure(fig, fullfile(cfg.output.dir, 'figures', 'ricker_wavelet_time_spectrum.png'));
close(fig);

write_time_csv(fullfile(cfg.output.dir, 'ricker_wavelet_time.csv'), t(:), wavelet(:), source(:));
write_spectrum_csv(fullfile(cfg.output.dir, 'ricker_wavelet_spectrum.csv'), freq(:), spectrum(:));
end

function write_time_csv(outFile, time, wavelet, source)
fid = fopen(outFile, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '时间_s,归一化振幅,物理震源振幅\n');
for k = 1:numel(time)
    fprintf(fid, '%.15g,%.15g,%.15g\n', time(k), wavelet(k), source(k));
end
end

function write_spectrum_csv(outFile, freq, spectrum)
fid = fopen(outFile, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '频率_Hz,归一化振幅谱\n');
for k = 1:numel(freq)
    fprintf(fid, '%.15g,%.15g\n', freq(k), spectrum(k));
end
end
