% fig2_wr_twist_sim.m
% Writhe and twist of the data and of the selected Cosserat-rod simulations
% (Fig. 2):
%   tw_plus_wr_10samples_V6.pdf   Tw vs Wr (data coloured by embryo, simulations
%                                 black), and Tw + Wr vs aligned time
%   wr_vs_time_sims_V6.pdf        Wr vs aligned time, data and simulations
% to P.out_fig2. The data-only version is Fig1/fig1_wr.m.
%
% The simulations start from an embryo's configuration at a given frame and
% run to the end of that embryo's recording; their saved states are spread
% uniformly over that frame range:
%   WT  : sample6_WT run, frames 54 .. end of sample 6, on sample 6's lag
%   Myo : sample1_myo run, frames 48 .. end of sample 3, on sample 3's lag
% (as in the published figure; the Myo run is placed on sample 3's timeline,
% see README).

clear; close all
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..')); P = config_paths();
if ~exist(P.out_fig2, 'dir'); mkdir(P.out_fig2); end

all_data = load(P.framesV6).all_data;
num_samples = numel(all_data);
time_lags = get_time_lags(P, all_data);

sims = struct( ...
    'file',     {'sample6_WT_T3.0_tw1.0_p10.0.mat', 'sample1_myo_T3.0_tw1.0_p10.0.mat'}, ...
    'f_start',  {54, 48}, ...
    'timeline', {6, 3});         % sample whose frames / lag the run is placed on

fig  = figure('Color', 'w', 'Position', [100, 100, 900, 420]);
fig2 = figure('Color', 'w');
for i = 1:num_samples
    d = all_data{i};
    [Wr, Tw] = compute_wr_tw(d.cent, d.Q);
    t_al = aligned_time((1:numel(Wr))', time_lags(i), 'frame', P);
    c = sample_color(i, num_samples);
    figure(fig)
    subplot(1,2,1); hold on; scatter(Wr, Tw, 'filled', 'MarkerFaceColor', c)
    subplot(1,2,2); hold on; plot(t_al, Tw + Wr, 'Color', c)
    figure(fig2); hold on
    plot(t_al, Wr, 'Color', c)
end

for k = 1:numel(sims)
    S = load(fullfile(P.sims, sims(k).file));
    [Wr, Tw] = compute_wr_tw(S.xsaves, S.Qsaves);
    j = sims(k).timeline;
    n_frames = size(all_data{j}.cent, 1);
    f = sims(k).f_start + (n_frames - sims(k).f_start) * linspace(0, 1, numel(Wr))';
    t_al = aligned_time(f, time_lags(j), 'frame', P);
    figure(fig)
    subplot(1,2,1); scatter(Wr, Tw, 'filled', 'MarkerFaceColor', 'k')
    subplot(1,2,2); plot(t_al, Tw + Wr, 'Color', 'k')
    figure(fig2)
    plot(t_al, Wr, 'Color', 'k', 'LineWidth', 2)
end

figure(fig)
subplot(1,2,1)
plot([-2, 2], -[-2, 2], 'Color', 'k')                       % Tw + Wr = 0, +-1
plot([-2, 2], -[-2, 2] + 1, 'Color', 'k', 'LineStyle', '--')
plot([-2, 2], -[-2, 2] - 1, 'Color', 'k', 'LineStyle', '--')
xlim([-1.5, 1.5]); ylim([-1.5, 1.5])
xlabel('Wr'); ylabel('Tw'); box on; axis square
subplot(1,2,2)
xlabel('Time [min]'); ylabel('Tw + Wr'); ylim([-1.5, 1.5])
box on; axis square
exportgraphics(fig, fullfile(P.out_fig2, 'tw_plus_wr_10samples_V6.pdf'), 'ContentType', 'vector')

figure(fig2)
xlabel('time (min)'); ylabel('Wr')
xlim([45, 120])
box on; axis square
fontsize(fig2, 9, 'points'); fontname(fig2, 'Helvetica')
exportgraphics(fig2, fullfile(P.out_fig2, 'wr_vs_time_sims_V6.pdf'), 'ContentType', 'vector')
