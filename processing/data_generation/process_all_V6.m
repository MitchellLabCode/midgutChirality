% process_all_V6.m
% Builds the framed-curve dataset all_framesInterpV6_10samples.mat from the
% TubULAR outputs (Sec. IV of the Supplement): smoothed centerline, myocyte
% stripes, material frame Q, Darboux vector kappa, spin omega, meshes, hoop
% radii and fold positions of every embryo. The structure of the file is
% described in Fig1/data_figshare_fig1/readme.txt.
%
% Inputs per sample (config_paths.m):
%   tubular_output/gridCoords_nU0100_nV0100/centerline_from_DVhoops_000step/  centerlines
%   tubular_output/cellTracking/muscleStripes/muscleNucleiStripes.mat          (extractNucleiStripes_noTracks.m)
%   tubular_output/gridCoords_nU0100_nV0100/features/fold_locations_sphi_nU0100_nV0100_avgpts.mat
%   P.meshes/meshes_sample<N>.mat                                              (tubi_to_meshrs.m)
% Output: P.framesV6

clear; close all
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..')); P = config_paths();

w = 0.3;          % Taubin / diffusive smoothing coefficient (space and time)
N_s = 20;         % smoothing iterations
n_samples = 100;  % points on the framed centerline
n_cut = 5;        % anterior centerline points discarded

all_data = cell(1, numel(P.samples));
for i = 1:numel(P.samples)
    dir_i = P.samples(i).tubular_output;
    fprintf('Processing sample #%d...\n', i);
    center_line_loc = fullfile(dir_i, 'gridCoords_nU0100_nV0100/centerline_from_DVhoops_000step');
    myocytes_loc    = fullfile(dir_i, 'cellTracking/muscleStripes/');
    folds_loc       = fullfile(dir_i, 'gridCoords_nU0100_nV0100/features/fold_locations_sphi_nU0100_nV0100_avgpts.mat');
    meshes_loc      = fullfile(P.meshes, ['meshes_sample' num2str(i)]);

    % smoothed centerline and myocyte stripes, then the material frame
    [x_raw, x] = process_sample_centerline(center_line_loc, n_samples + n_cut, w, N_s);
    M = load(meshes_loc);
    [myocytes_raw, myocytes] = load_myocytes_stripes(myocytes_loc, n_samples, w, N_s);
    [cent, Q, kappa, omega] = makeFramedCurve(x, myocytes, M.meshes, n_cut, w, N_s);

    F = load(folds_loc);
    all_data{i} = struct('dir', string(dir_i), 'x_raw', x_raw, 'x', x, ...
        'myocytes_raw', myocytes_raw, 'myocytes', myocytes, 'cent', cent, 'Q', Q, ...
        'kappa', kappa, 'omega', omega, 'meshes', {M.meshes}, 'radii', M.radii, ...
        'ssmax', F.ssmax, 'folds', F.folds, 'ssfold', F.ssfold, ...
        'fold_onset', F.fold_onset, 'timeInterval', F.timeInterval);
end
save(P.framesV6, 'all_data', '-v7.3')
fprintf('wrote %s\n', P.framesV6);
