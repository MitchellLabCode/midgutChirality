% export_figshare_fig2.m
% Figshare data for Fig. 2, in P.figshare_fig2, from all_framesInterpV6_10samples.mat:
%
%   fig2_kymographs.h5          per embryo: radius R(s,t), contact radius R_c(s,t),
%                               edge growth rate, and fold / A-P split positions
%   fig2_curvature_measures.csv E_b1, E_b2, E_b, E_tw (Sec. IV, Eqs. S16-S18), per unit
%                               length and integrated, per region
%   fig2_rotation_rates.csv     average rotation speed omega_b, omega_tw (Eq. S15), and the
%                               d1 / d2 components omega_1, omega_2, per region
%   fig2_curvature_change.csv   average curvature change mu_b, mu_tw (Eq. S19), and mu_1, mu_2
%   readme.txt                  copy of figshare_readme.txt
%
% The measures are those plotted by fig2_deformation_measures.m
% (compute_deformation_measures.m), the radius by fig2_radius.m, the growth
% rate by fig2_growthrate.m and the contact radius by
% bounding_contact/generate_bounding_contact.m (contact_radius.m).
% Regions: anterior, posterior, all_midgut.

clear; close all
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..')); P = config_paths();

out_dir = P.figshare_fig2;
h5file  = fullfile(out_dir, 'fig2_kymographs.h5');
dt = P.dt;

all_data = load(P.framesV6).all_data;
num_samples = numel(all_data);
time_lags = get_time_lags(P, all_data);

if ~exist(out_dir, 'dir'); mkdir(out_dir); end
if exist(h5file, 'file'); delete(h5file); end
regions = {'anterior', 'posterior', 'all_midgut'};
Tcurv = {}; Trot = {}; Tchg = {};

for idx = 1:num_samples
    d = all_data{idx};
    x = d.cent;
    num_t = size(d.Q, 1);
    is_wt = P.samples(idx).is_WT;
    lag = time_lags(idx);

    %%% deformation measures
    D = compute_deformation_measures(d, lag, P);
    n = numel(D.E.t);
    for r = 1:3
        Tcurv{end+1} = table(D.E.t, rep(idx, n), rep(is_wt, n), rep(string(regions{r}), n), ...
            D.E.E_b1(:,r), D.E.E_b2(:,r), D.E.E_b(:,r), D.E.E_tw(:,r), ...
            D.E.E_b1_total(:,r), D.E.E_b2_total(:,r), D.E.E_b_total(:,r), D.E.E_tw_total(:,r), ...
            'VariableNames', {'time_min','sample_id','is_WT','region', ...
            'E_b1','E_b2','E_b','E_tw','E_b1_total','E_b2_total','E_b_total','E_tw_total'}); %#ok<SAGROW>
    end
    n = numel(D.omega.t);
    for r = 1:3
        Trot{end+1} = table(D.omega.t, rep(idx, n), rep(is_wt, n), rep(string(regions{r}), n), ...
            D.omega.omega_b(:,r), D.omega.omega_tw(:,r), D.omega.omega_1(:,r), D.omega.omega_2(:,r), ...
            'VariableNames', {'time_min','sample_id','is_WT','region','omega_b','omega_tw','omega_1','omega_2'}); %#ok<SAGROW>
    end
    n = numel(D.mu.t);
    for r = 1:3
        Tchg{end+1} = table(D.mu.t, rep(idx, n), rep(is_wt, n), rep(string(regions{r}), n), ...
            D.mu.mu_b(:,r), D.mu.mu_tw(:,r), D.mu.mu_1(:,r), D.mu.mu_2(:,r), ...
            'VariableNames', {'time_min','sample_id','is_WT','region','mu_b','mu_tw','mu_1','mu_2'}); %#ok<SAGROW>
    end

    %%% kymographs
    nm = min(num_t, numel(d.meshes));
    hs = hoop_arclength(d); hs = hs(1:nm, :);      % s/L of each hoop, per frame
    Rc = contact_radius(d);
    ell = vecnorm(diff(x, 1, 2), 2, 3);            % (num_t, N-1) edge lengths
    growth = diff(ell, 1, 1) / dt;                 % (num_t-1, N-1) [um/min]
    fs_cent = fold_arclength(d, 'cent');
    nf = min(num_t, size(d.folds, 1));
    fold_hoop = nan(num_t, 3); fold_s = nan(num_t, 3);
    for k = 1:3
        ok = find((1:nf)' > d.fold_onset(k));
        fold_hoop(ok, k) = double(d.folds(ok, k));
        fold_s(ok, k) = fs_cent(ok, k);
    end
    [~, idx_split] = max(x(:,:,1), [], 2);
    split_s = (idx_split - 1) / (size(x,2) - 1);

    g = sprintf('/sample_%02d', idx);
    w5 = @(name, val) write_ds(h5file, [g '/' name], val);
    w5('R', d.radii);                              w5('R_s_over_L', hs);
    w5('R_time_min', aligned_time((1:num_t)', lag, 'frame', P));
    w5('R_contact', Rc);                           w5('cent_s_over_L', linspace(0, 1, size(x,2))');
    w5('edge_growth_rate', growth);                w5('edge_s_over_L', ((1:size(ell,2))' - 0.5) / size(ell,2));
    w5('edge_growth_time_min', aligned_time((1:num_t-1)', lag, 'interval', P));
    w5('fold_hoop_index', fold_hoop);              w5('fold_s_over_L', fold_s);
    w5('fold_onset_frame', d.fold_onset(:));       w5('max_x_s_over_L', split_s);
    h5writeatt(h5file, g, 'sample_id', int32(idx));
    h5writeatt(h5file, g, 'is_WT', uint8(is_wt));
    h5writeatt(h5file, g, 'genotype', P.samples(idx).genotype);
    h5writeatt(h5file, g, 'time_lag_frames', lag);
    h5writeatt(h5file, g, 'time_lag_min', dt * lag);
    h5writeatt(h5file, g, 'n_frames', int32(num_t));
    h5writeatt(h5file, g, 'source_dir', char(d.dir));
    fprintf('sample %2d done\n', idx);
end

h5writeatt(h5file, '/', 'description', ['Fig. 2 kymographs per embryo: mesh-hoop radius R, contact radius ' ...
    'R_contact, edge growth rate, fold and A/P-split positions. See readme.txt.']);
h5writeatt(h5file, '/', 'data_file', 'all_framesInterpV6_10samples.mat');
h5writeatt(h5file, '/', 'frame_interval_min', dt);
h5writeatt(h5file, '/', 'created', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));

writetable(vertcat(Tcurv{:}), fullfile(out_dir, 'fig2_curvature_measures.csv'));
writetable(vertcat(Trot{:}),  fullfile(out_dir, 'fig2_rotation_rates.csv'));
writetable(vertcat(Tchg{:}),  fullfile(out_dir, 'fig2_curvature_change.csv'));
copyfile(fullfile(here, 'figshare_readme.txt'), fullfile(out_dir, 'readme.txt'));
fprintf('wrote %s and 3 csv files\n', h5file);


function c = rep(v, n)
c = repmat(v, n, 1);
end

function write_ds(f, name, val)
% vectors are written as 1-D datasets (shape (n,) in h5py), matrices as 2-D
if isvector(val), sz = numel(val); val = val(:); else, sz = size(val); end
h5create(f, name, sz, 'Datatype', 'double');
h5write(f, name, double(val));
end
