% prep_interp_data.m
% Inputs of the canalization simulations (cosserat/canalization/): shapes
% interpolated between a WT and a byn>Myo1C embryo (Fig. S15).
%
% For 50 values of tau from -0.1 to 1.1 (tau = 0: WT, tau = 1: byn>Myo1C) and
% every frame common to both embryos from aligned frame t0 on:
%   x_i          : centerlines, SRVinterpolation(x1, x2, tau) + (1-tau) x1 + tau x2
%   Q_i          : frames: the byn>Myo1C frame is minimally rotated onto the WT
%                  tangent, the WT frame is rotated about that tangent by tau
%                  times the angle between them, and the result is minimally
%                  rotated onto the tangent of x_i (FrameInterpolation.m)
%   circum_radii : (1-tau) x the WT contact radius + tau x the byn>Myo1C one
%   x0           : the same interpolation of the two start shapes
%                  (cosserat/prod/x0_*.mat)
% Walls n, d (n . x < d): the convex hull of all x_i (every tau and frame) and
% all x0, the same for every tau.
%
% Reads  P.bounding_contact/refIV6_sample<N>_walls3_cball.mat (generate_bounding_contact.m)
% Writes P.canalization/SRVbox3f_s<WT>s<Myo>/itp_<k>_tau<tau>.mat with x_i,
%        Q_i, circum_radii, n, d, time_lags (all samples), tau and x0.

clear; close all
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'), here, fullfile(here, '..', 'bounding_contact'));
P = config_paths();

% samples to interpolate, and their start shapes
idxWT  = 6;
idxMyo = 3;
x0_WT  = fullfile(here, '..', '..', 'cosserat', 'prod', 'x0_s6V5_t54_ant.mat');
x0_Myo = fullfile(here, '..', '..', 'cosserat', 'prod', 'x0_s3V5_t45_f3.mat');

t0 = 19;                                % first aligned frame (frame + lag)
taus = linspace(-0.1, 1.1, 50);
show_plots = false;                     % 3D view of the points spanning the walls

out_dir = fullfile(P.canalization, sprintf('SRVbox3f_s%ds%d', idxWT, idxMyo));

%%% reference data
ref = @(k) fullfile(P.bounding_contact, sprintf('refIV6_sample%d_walls3_cball.mat', k));
time_lags = zeros(1, numel(P.samples));
for k = 1:numel(P.samples)
    time_lags(k) = load(ref(k), 'tlag').tlag;
end
dat_WT  = load(ref(idxWT));
dat_Myo = load(ref(idxMyo));
x_wt0  = load(x0_WT, 'x0').x0;
x_myo0 = load(x0_Myo, 'x0').x0;

% frames common to both embryos, from t0 to the end of the shorter series
time_arr_wt  = (1:size(dat_WT.x, 1))  + time_lags(idxWT);
time_arr_myo = (1:size(dat_Myo.x, 1)) + time_lags(idxMyo);
tf = min(time_arr_myo(end), time_arr_wt(end));
tpsWT  = find(time_arr_wt  >= t0 & time_arr_wt  <= tf);
tpsMyo = find(time_arr_myo >= t0 & time_arr_myo <= tf);
num_tp = numel(tpsMyo);

x_wt  = dat_WT.x;   Q_wt  = dat_WT.Q;   cr_wt  = dat_WT.circum_radii(tpsWT,:);
x_myo = dat_Myo.x;  Q_myo = dat_Myo.Q;  cr_myo = dat_Myo.circum_radii(tpsMyo,:);

%%% interpolate at every tau (kept, so the walls can enclose all of them)
num_tau = numel(taus);
x_all = cell(1, num_tau); Q_all = cell(1, num_tau);
cr_all = cell(1, num_tau); x0_all = cell(1, num_tau);
for k = 1:num_tau
    tau = taus(k);
    x_i = zeros(num_tp, size(x_wt, 2), 3);
    Q_i = zeros(num_tp, size(Q_wt, 2), 3, 3);
    for t = 1:num_tp
        x1 = squeeze(x_wt(tpsWT(t),:,:));
        x2 = squeeze(x_myo(tpsMyo(t),:,:));
        x_i(t,:,:) = SRVinterpolation(x1, x2, tau) + (1-tau)*x1 + tau*x2;
        Q_i(t,:,:,:) = FrameInterpolation(squeeze(x_i(t,:,:)), ...
            squeeze(Q_wt(tpsWT(t),:,:,:)), squeeze(Q_myo(tpsMyo(t),:,:,:)), tau);
    end
    x_all{k} = x_i;
    Q_all{k} = Q_i;
    cr_all{k} = (1-tau)*cr_wt + tau*cr_myo;
    x0_all{k} = SRVinterpolation(x_wt0, x_myo0, tau) + (1-tau)*x_wt0 + tau*x_myo0;
end

%%% walls: convex hull of every interpolated centerline and start shape
pts = [cell2mat(cellfun(@(x) reshape(x, [], 3), x_all(:), 'UniformOutput', false));
       vertcat(x0_all{:})];
[n, d] = boundingplanes(pts(:,1), pts(:,2), pts(:,3));
% convhull 'Simplify' can drop near-coplanar vertices, leaving a few points
% slightly outside; push each plane out so that every point is inside
d = max(d, max(pts * n.', [], 1).');

if show_plots
    figure('Color', 'w')
    scatter3(pts(:,1), pts(:,2), pts(:,3), 1)
    daspect([1,1,1])
end

%%% save
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
for k = 1:num_tau
    tau = taus(k);
    x_i = x_all{k}; Q_i = Q_all{k}; circum_radii = cr_all{k}; x0 = x0_all{k};
    save(fullfile(out_dir, ['itp_' num2str(k, '%02d') '_tau' num2str(tau) '.mat']), ...
        'x_i', 'Q_i', 'circum_radii', 'n', 'd', 'time_lags', 'tau', 'x0');
end
fprintf('wrote %d files to %s\n', num_tau, out_dir);
