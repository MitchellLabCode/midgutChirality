% fig2_deformation_measures.m
% Time-aligned deformation measures of the framed centerline (Fig. 2 and
% Supplement Sec. IV.4), per region (anterior, posterior, whole gut), for all
% 10 embryos: thin lines per embryo (WT solid, byn>Myo1C dashed), thick
% lines the binned mean. All measures come from compute_deformation_measures.m
% (median filtering, Voronoi / edge-length weights, A/P split at the max-x
% point); this script only selects and plots them. Writes to P.out_fig2:
%
%   timealigned_kappaEnergy_10samples_split.pdf        E_tw, E_b: mean, integral, total angle
%   timealigned_bendEnergy_d1d2_10samples_split.pdf    E_b1, E_b2: mean, integral, total angle
%   timealigned_intEnergy_d1d2d3_10samples_split.pdf   integrals of E_b1, E_b2, E_tw
%   timealigned_darboux_10samples_split.pdf            mu_tw, mu_b
%   timealigned_darboux_d1d2_twist_10samples_split.pdf mu_1, mu_2, mu_tw
%   timealigned_twist_bend_10samples_voro_split.pdf    omega_tw, omega_b
%   timealigned_spin_d1d2_twist_10samples_split.pdf    omega_1, omega_2, omega_tw
%   timealigned_dtkappa_omega_10samples_split.pdf      mu_tw, mu_b and omega_tw, omega_b
%
% Notation (Supplement Sec. IV): kappa_i, omega_i are the material components
% of the Darboux vector and the frame spin; E_* = kappa_*^2/2; mu = |d_t kappa|.

clear; close all
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..')); P = config_paths();
if ~exist(P.out_fig2, 'dir'); mkdir(P.out_fig2); end

all_data = load(P.framesV6).all_data;
num_samples = numel(all_data);
time_lags = get_time_lags(P, all_data);

S = struct('is_WT', {}, 'E', {}, 'omega', {}, 'mu', {});
for i = 1:num_samples
    D = compute_deformation_measures(all_data{i}, time_lags(i), P);
    S(i).is_WT = P.samples(i).is_WT;
    S(i).E = D.E; S(i).omega = D.omega; S(i).mu = D.mu;
end

gem = orderedcolors('gem');
c1 = gem(1,:); c2 = gem(2,:); c3 = gem(3,:);
red = [1 0 0]; blk = [0 0 0];
ylab_E     = 'Mean \kappa^2/2 [\mum^{-2}]';
ylab_Eint  = '\int\kappa^2/2 ds [\mum^{-1}]';
ylab_Theta = '\int|\kappa| ds [rad]';
ylab_mu    = 'Mean |\partial_t\kappa| [rad/(\mum\cdotmin)]';
ylab_om    = 'Mean rotation speed [rad/min]';
row = @(src, names, colors, labels, ylab) struct('src', src, 'names', {names}, ...
    'colors', {colors}, 'labels', {labels}, 'ylabel', ylab, 'ylim_from', 'all');
tw_b = {'twist', 'bend'};

figs = struct('name', {}, 'rows', {});
figs(end+1) = struct('name', 'timealigned_kappaEnergy_10samples_split', 'rows', [ ...
    row('E', {'E_tw', 'E_b'}, {red, blk}, tw_b, ylab_E), ...
    row('E', {'E_tw_total', 'E_b_total'}, {red, blk}, tw_b, ylab_Eint), ...
    row('E', {'Theta_tw', 'Theta_b'}, {red, blk}, tw_b, ylab_Theta)]);
figs(end+1) = struct('name', 'timealigned_bendEnergy_d1d2_10samples_split', 'rows', [ ...
    row('E', {'E_b1', 'E_b2'}, {c1, c2}, {'\kappa_1', '\kappa_2'}, ylab_E), ...
    row('E', {'E_b1_total', 'E_b2_total'}, {c1, c2}, {'\kappa_1', '\kappa_2'}, ylab_Eint), ...
    row('E', {'Theta_b1', 'Theta_b2'}, {c1, c2}, {'\kappa_1', '\kappa_2'}, ylab_Theta)]);
r = row('E', {'E_b1_total', 'E_b2_total', 'E_tw_total'}, {c1, c2, c3}, ...
    {'bend d_1', 'bend d_2', 'twist d_3'}, '\int\kappa_i^2/2 ds [\mum^{-1}]');
r.ylim_from = 'whole';      % whole-gut integral = anterior + posterior: sets the scale
figs(end+1) = struct('name', 'timealigned_intEnergy_d1d2d3_10samples_split', 'rows', r);
figs(end+1) = struct('name', 'timealigned_darboux_10samples_split', 'rows', ...
    row('mu', {'mu_tw', 'mu_b'}, {red, blk}, tw_b, ylab_mu));
figs(end+1) = struct('name', 'timealigned_darboux_d1d2_twist_10samples_split', 'rows', ...
    row('mu', {'mu_1', 'mu_2', 'mu_tw'}, {c1, c2, c3}, {'\kappa_1', '\kappa_2', 'twist'}, ylab_mu));
figs(end+1) = struct('name', 'timealigned_twist_bend_10samples_voro_split', 'rows', ...
    row('omega', {'omega_tw', 'omega_b'}, {red, blk}, tw_b, ylab_om));
figs(end+1) = struct('name', 'timealigned_spin_d1d2_twist_10samples_split', 'rows', ...
    row('omega', {'omega_1', 'omega_2', 'omega_tw'}, {c1, c2, c3}, {'\omega_1', '\omega_2', 'twist'}, ylab_om));
figs(end+1) = struct('name', 'timealigned_dtkappa_omega_10samples_split', 'rows', [ ...
    row('mu', {'mu_tw', 'mu_b'}, {red, blk}, tw_b, ylab_mu), ...
    row('omega', {'omega_tw', 'omega_b'}, {red, blk}, tw_b, 'Mean rotation speed |\omega| [rad/min]')]);

for k = 1:numel(figs)
    fprintf('%s\n', figs(k).name);
    fig = plot_timealigned(S, figs(k).rows);
    exportgraphics(fig, fullfile(P.out_fig2, [figs(k).name '.pdf']), 'ContentType', 'vector')
    close(fig)
end
