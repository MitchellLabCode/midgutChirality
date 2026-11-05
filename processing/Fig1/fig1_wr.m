% fig1_wr.m
% Writhe Wr of the centerline and twist Tw of the frame (d1 about the
% centerline) for all 10 embryos, on the aligned time axis (Fig. 1).
%
%   tw_wr_10samples_V6.pdf           Wr(t) and Tw(t)
%   tw_plus_wr_vs_time_10samples_V6.pdf   Tw + Wr vs time
%   tw_plus_wr_10samples_data_V6.pdf      Tw vs Wr, and Tw + Wr vs time (data only;
%                                         the version with simulations is Fig2/fig2_wr_twist_sim.m)
%   wr_tw_10samples_V6.csv           aligned_time_min, sample_id, is_WT, Wr, Tw
% to P.out_fig1. Wr and Tw are computed by compute_writhe.m and compute_twist.m.

clear; close all
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..')); P = config_paths();
if ~exist(P.out_fig1, 'dir'); mkdir(P.out_fig1); end

all_data = load(P.framesV6).all_data;
num_samples = numel(all_data);
time_lags = get_time_lags(P, all_data);

fig1 = figure('Color', 'w');     % Wr, Tw vs time
fig2 = figure('Color', 'w');     % Tw + Wr vs time
fig3 = figure('Color', 'w', 'Position', [100, 100, 900, 420]);   % Tw vs Wr, Tw + Wr vs time
csv_rows = cell(1, num_samples);
for idx_sample = 1:num_samples
    d = all_data{idx_sample};
    [Wr, Tw] = compute_wr_tw(d.cent, d.Q);
    num_t = numel(Wr);
    t_al = aligned_time((1:num_t)', time_lags(idx_sample), 'frame', P);
    c = sample_color(idx_sample, num_samples);

    csv_rows{idx_sample} = table(t_al, repmat(idx_sample, num_t, 1), ...
        repmat(P.samples(idx_sample).is_WT, num_t, 1), Wr, Tw, ...
        'VariableNames', {'aligned_time_min', 'sample_id', 'is_WT', 'Wr', 'Tw'});

    figure(fig1)
    subplot(2,1,1); hold on; plot(t_al, Wr, 'Color', c)
    subplot(2,1,2); hold on; plot(t_al, Tw, 'Color', c)
    figure(fig2); hold on
    plot(t_al, Tw + Wr, 'Color', c)
    figure(fig3)
    subplot(1,2,1); hold on; scatter(Wr, Tw, 'filled', 'MarkerFaceColor', c)
    subplot(1,2,2); hold on; plot(t_al, Tw + Wr, 'Color', c)
end

figure(fig1)
subplot(2,1,1); ylabel('Wr'); ylim([-1.5, 1.5]); xlim([-50, 150]); box on
subplot(2,1,2); xlabel('Time [min]'); ylabel('Tw'); ylim([-1.5, 1.5]); xlim([-50, 150]); box on

figure(fig2)
xlabel('Time [min]'); ylabel('Tw + Wr'); ylim([-1.5, 1.5])
yline(0, 'Color', 'k', 'LineStyle', '--')
box on; axis square

figure(fig3)
subplot(1,2,1)
plot([-2, 2], -[-2, 2], 'Color', 'k')                       % Tw + Wr = 0, +-1
plot([-2, 2], -[-2, 2] + 1, 'Color', 'k', 'LineStyle', '--')
plot([-2, 2], -[-2, 2] - 1, 'Color', 'k', 'LineStyle', '--')
xlim([-1.5, 1.5]); ylim([-1.5, 1.5])
xlabel('Wr'); ylabel('Tw'); box on; axis square
subplot(1,2,2)
xlabel('Time [min]'); ylabel('Tw + Wr'); ylim([-1.5, 1.5]); xlim([-50, 150])
box on; axis square

exportgraphics(fig1, fullfile(P.out_fig1, 'tw_wr_10samples_V6.pdf'), 'ContentType', 'vector')
exportgraphics(fig2, fullfile(P.out_fig1, 'tw_plus_wr_vs_time_10samples_V6.pdf'), 'ContentType', 'vector')
exportgraphics(fig3, fullfile(P.out_fig1, 'tw_plus_wr_10samples_data_V6.pdf'), 'ContentType', 'vector')
T = vertcat(csv_rows{:});
writetable(T, fullfile(P.out_fig1, 'wr_tw_10samples_V6.csv'))
fprintf('wrote Wr/Tw figures and wr_tw_10samples_V6.csv (%d rows) to %s\n', height(T), P.out_fig1);
