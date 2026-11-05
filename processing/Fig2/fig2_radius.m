% fig2_radius.m
% Gut radius (Fig. 2). Top: mean hoop radius vs aligned time for all
% embryos. Bottom: radius kymograph R(s/L, t) of the reference embryo
% (sample 4), Gaussian-filtered (sigma = 1), on the arclength of each hoop
% centroid, with the folds (grey) and the max-x point (dashed).
% Writes fig2_radius.pdf to P.out_fig2.

clear; close all
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..')); P = config_paths();
if ~exist(P.out_fig2, 'dir'); mkdir(P.out_fig2); end

all_data = load(P.framesV6).all_data;
num_samples = numel(all_data);
time_lags = get_time_lags(P, all_data);
idx_ref = 4;

fig = figure('Color', 'w');
tiledlayout(4, 3);
nexttile([1, 3]); hold on
for i = 1:num_samples
    mean_rad = mean(all_data{i}.radii, 2);
    plot(aligned_time((1:numel(mean_rad))', time_lags(i), 'frame', P), mean_rad, 'Color', sample_color(i, num_samples))
end
ylabel('mean radius [\mum]'); box on

d = all_data{idx_ref};
radius = d.radii;
num_t = size(radius, 1);
t_al = aligned_time((1:num_t)', time_lags(idx_ref), 'frame', P);
xlim([t_al(1) t_al(end)])

nexttile(4, [3, 3])
colormap(bone)
radius = imgaussfilt(radius, 1);
% hoops are not evenly spaced in arclength: put them on s/L of the hoop centroids
hs = hoop_arclength(d);
s_grid = linspace(0, 1, size(radius, 2));
radius_s = zeros(size(radius));
for t = 1:num_t
    radius_s(t,:) = interp1(hs(t,:), radius(t,:), s_grid, 'linear');
end
imagesc('XData', t_al, 'YData', s_grid, 'CData', radius_s')
hold on
folds = d.folds;
for k = 1:3, folds(1:d.fold_onset(k), k) = -1; end
for k = 1:3     % folds are hoop indices: place them at that hoop's s/L
    valid = find(folds(:,k) > 0 & (1:size(folds,1))' <= size(hs,1));
    plot(t_al(valid), hs(sub2ind(size(hs), valid, double(folds(valid,k)))), 'Color', 0.6*[1,1,1], 'LineWidth', 1)
end
% max-x point of cent, as s/L of the full centerline x
n_cut = size(d.x, 2) - size(d.cent, 2);
nt_c = size(d.cent, 1);
split_s = zeros(nt_c, 1);
for t = 1:nt_c
    [~, imax] = max(d.cent(t,:,1));
    ss = [0; cumsum(vecnorm(diff(squeeze(d.x(t,:,:))), 2, 2))];
    split_s(t) = ss(imax + n_cut) / ss(end);
end
plot(t_al(1:nt_c), split_s, 'Color', 0.6*[1,1,1], 'LineStyle', '--')
ylim([0, 1]); xlim([t_al(1) t_al(end)])
xlabel('Time [min]'); ylabel('s/L')
cb = colorbar; cb.Label.String = 'R [\mum]';
exportgraphics(fig, fullfile(P.out_fig2, 'fig2_radius.pdf'), 'ContentType', 'vector')
