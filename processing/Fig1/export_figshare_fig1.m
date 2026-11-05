% export_figshare_fig1.m
% Figshare data for Fig. 1, in P.figshare_fig1:
%   chiral_kymographs.h5     the chirality kymographs q(s, t) of all 10 embryos,
%                            computed exactly as in fig1_chirality_10samples.m
%   wr_tw_10samples_V6.csv   writhe and twist (written by fig1_wr.m; run it first)
%   readme.txt               copy of figshare_readme.txt
%
% HDF5 layout (one group per sample):
%   /                         attrs: description, data_file, w, n_skip_start, n_skip_end, ...
%   /sample_NN                attrs: sample_id, is_WT, genotype, time_lag_frames, time_lag_min,
%                                    n_frames, source_dir
%   /sample_NN/q              (n_s x n_frames)  chirality q(s, t)
%   /sample_NN/s_over_L       (n_s)             arclength fraction of each row of q
%   /sample_NN/time_min       (n_frames)        aligned time of each column [min]
%   /sample_NN/fold_s_over_L  (n_frames x 3)    fold positions (s/L), NaN before onset
%   /sample_NN/max_x_s_over_L (n_frames)        s/L of the max-x point (A/P split)

clear; close all
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..')); P = config_paths();

out_dir = P.figshare_fig1;
h5file  = fullfile(out_dir, 'chiral_kymographs.h5');
w = 5; n_skip_start = 9; n_skip_end = 1;

all_data = load(P.framesV6).all_data;
num_samples = numel(all_data);
time_lags = get_time_lags(P, all_data);

if ~exist(out_dir, 'dir'); mkdir(out_dir); end
if exist(h5file, 'file'); delete(h5file); end

for idx = 1:num_samples
    data = all_data{idx};
    C = arclength_resample(data.x);
    num_t = size(C, 1);
    [q, s_kymo] = coiling_kymograph(C, w, n_skip_start, n_skip_end);

    split_s = zeros(num_t, 1);
    for t = 1:num_t
        ce = squeeze(C(t,:,:)); ss = [0; cumsum(vecnorm(diff(ce), 2, 2))];
        [~, i] = max(ce(:,1)); split_s(t) = ss(i) / ss(end);
    end
    fold_s = nan(num_t, 3);
    fs = fold_arclength(data, 'x');
    nf = min(num_t, size(fs, 1));
    for k = 1:3
        valid = find((1:nf)' > data.fold_onset(k) & data.folds(1:nf,k) > 0);
        fold_s(valid, k) = fs(valid, k);
    end
    tps = aligned_time((1:num_t)', time_lags(idx), 'frame0', P);

    g = sprintf('/sample_%02d', idx);
    h5create(h5file, [g '/q'], size(q));                     h5write(h5file, [g '/q'], q);
    h5create(h5file, [g '/s_over_L'], numel(s_kymo));        h5write(h5file, [g '/s_over_L'], s_kymo(:));
    h5create(h5file, [g '/time_min'], num_t);                h5write(h5file, [g '/time_min'], tps(:));
    h5create(h5file, [g '/fold_s_over_L'], size(fold_s));    h5write(h5file, [g '/fold_s_over_L'], fold_s);
    h5create(h5file, [g '/max_x_s_over_L'], num_t);          h5write(h5file, [g '/max_x_s_over_L'], split_s);

    is_wt = P.samples(idx).is_WT;
    if is_wt, geno = 'WT (HandGFP;byn-GAL4,klar > UAS-mCh-CAAX,H2Av-iRFP)';
    else,     geno = 'byn>Myo1C (HandGFP;byn-GAL4,klar > UAS-Myo1C-RFP,H-iRFP)'; end
    h5writeatt(h5file, g, 'sample_id', int32(idx));
    h5writeatt(h5file, g, 'is_WT', uint8(is_wt));
    h5writeatt(h5file, g, 'genotype', geno);
    h5writeatt(h5file, g, 'time_lag_frames', time_lags(idx));
    h5writeatt(h5file, g, 'time_lag_min', P.dt * time_lags(idx));
    h5writeatt(h5file, g, 'n_frames', int32(num_t));
    h5writeatt(h5file, g, 'source_dir', char(data.dir));
    fprintf('sample %2d (%s): q %s, lag %+g frames, t %g..%g min\n', idx, ...
        P.samples(idx).genotype, mat2str(size(q)), time_lags(idx), tps(1), tps(end));
end

if strcmp(P.time_convention, 'unified')
    tconv = 'time_min = 2*frame + 2*time_lag_frames';
else
    tconv = 'time_min = 2*(frame-1) + 2*time_lag_frames';
end
h5writeatt(h5file, '/', 'description', ['Quasi-local chirality measure q(s, t) of the midgut ' ...
    'centerline: signed solid angle swept by the unit tangent over a sliding window of 2w+1 tangents ' ...
    '(qmeasure.m). One group per embryo. See readme.txt.']);
h5writeatt(h5file, '/', 'data_file', 'all_framesInterpV6_10samples.mat');
h5writeatt(h5file, '/', 'window_half_width_w', int32(w));
h5writeatt(h5file, '/', 'n_skip_start', int32(n_skip_start));
h5writeatt(h5file, '/', 'n_skip_end', int32(n_skip_end));
h5writeatt(h5file, '/', 'frame_interval_min', P.dt);
h5writeatt(h5file, '/', 'time_convention', [tconv '; t = 0 at the onset ' ...
    'of the 2nd fold of the reference WT embryo (sample 4)']);
h5writeatt(h5file, '/', 'created', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
fprintf('wrote %s\n', h5file);

csv = fullfile(P.out_fig1, 'wr_tw_10samples_V6.csv');
if exist(csv, 'file')
    copyfile(csv, out_dir);
else
    warning('%s not found: run fig1_wr.m first', csv);
end
copyfile(fullfile(here, 'figshare_readme.txt'), fullfile(out_dir, 'readme.txt'));
