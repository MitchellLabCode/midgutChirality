% fig2_growthrate.m
% Edge growth kymographs (Fig. 2): for every embryo, the change of each
% centerline edge length between consecutive frames, l(t+1) - l(t) [um per
% frame], on s/L of cent and the aligned time of each interval, with the
% folds (grey) and the max-x point (dashed). One panel per embryo, titled
% '<sample> . <genotype>'. Writes fig2_growthrate.pdf to P.out_fig2.

clear; close all
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..')); P = config_paths();
if ~exist(P.out_fig2, 'dir'); mkdir(P.out_fig2); end

all_data = load(P.framesV6).all_data;
num_samples = numel(all_data);
time_lags = get_time_lags(P, all_data);

fig = figure('Color', 'w', 'Position', [100, 100, 1500, 700]);
tiledlayout(2, ceil(num_samples/2));
for i = 1:num_samples
    d = all_data{i};
    x_all = d.cent;
    [num_t, num_s, ~] = size(x_all);
    lengths = vecnorm(diff(x_all, 1, 2), 2, 3);          % (n_t x N-1) edge lengths
    growth = diff(lengths, 1, 1);                        % (n_t-1 x N-1)
    t_int = aligned_time((1:num_t-1)', time_lags(i), 'interval', P);
    t_frm = aligned_time((1:num_t)', time_lags(i), 'frame', P);

    nexttile
    imagesc('XData', [0, 1], 'YData', t_int, 'CData', growth)
    hold on
    colormap(blues())
    folds = d.folds;
    for k = 1:3, folds(1:d.fold_onset(k), k) = -1; end
    fold_s = fold_arclength(d, 'cent');
    for k = 1:3
        valid = folds(:,k) > 0;
        plot(fold_s(valid,k), t_frm(valid), 'Color', 0.6*[1,1,1], 'LineWidth', 1)
    end
    [~, idx_split] = max(x_all(:,:,1), [], 2);
    plot((idx_split - 1)/(num_s - 1), t_frm, 'Color', 0.6*[1,1,1], 'LineStyle', '--')
    xlim([0, 1]); ylim([t_frm(1), t_frm(end)])
    colorbar; clim([-0.1, 0.6])
    box on
    set(gca, 'YDir', 'reverse')
    title(sprintf('%d \\cdot %s', i, P.samples(i).genotype))
    if i == 1 || i == ceil(num_samples/2) + 1, ylabel('Time [min]'); end
    if i > ceil(num_samples/2), xlabel('s/L'); end
end
exportgraphics(fig, fullfile(P.out_fig2, 'fig2_growthrate.pdf'), 'ContentType', 'vector')
