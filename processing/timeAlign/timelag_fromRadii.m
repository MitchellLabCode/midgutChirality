% timelag_fromRadii.m
% Time alignment of all embryos from their radius kymographs R(s, t).
%
% 1. For every embryo, the mean hoop radius of each cut mesh is read from
%    TubULAR and interpolated onto a fold-anchored arclength coordinate
%    (interpolate_spatially_between_folds.m), so that kymographs of
%    different embryos are comparable position by position.
% 2. Within each genotype, every embryo is cross-correlated in time with a
%    reference embryo (compute_radius_timestamps.m ->
%    compute_kymoshifted_temporal_xcorrelation.m): rows are z-scored, the
%    Pearson correlation of the overlapping rows is computed for every integer
%    frame shift (|shift| <= 30), and the optimum of the Savitzky-Golay
%    smoothed correlation curve is the shift.
% 3. The two reference embryos are cross-correlated to give one cross-genotype
%    shift (compute_cross_genotype_timing.m).
%
% Writes P.timeshifts (timeshifts_fromRadiusKymos_10samples.mat), read by
% get_time_lags.m in every plotting script, and the cross-correlation figures
% (to output/ next to it).
% Requires TubULAR, and tubi.identifyFolds to have been run for every embryo.
%
% Internal numbering: datasets 1-5 are the WT embryos (samples 2, 4, 6, 8, 10),
% datasets 6-10 the byn>Myo1C embryos (samples 1, 3, 5, 7, 9); the reference
% embryos are dataset 2 (WT #2 = sample 4) and dataset 6 (byn>Myo1C #1 = sample 1).

clear; close all
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'), here, fullfile(here, '..', 'data_generation')); P = config_paths();
outdir = fullfile(fileparts(P.timeshifts), 'output');   % cross-correlation figures
if ~exist(outdir, 'dir'); mkdir(outdir); end
overwrite = true;

wt_samples  = find([P.samples.is_WT]);
myo_samples = find(~[P.samples.is_WT]);
dset_samples = [wt_samples, myo_samples];
genotypes = [repmat({'wt'}, 1, numel(wt_samples)), repmat({'bynMyo1c'}, 1, numel(myo_samples))];
dataDirs = {P.samples(dset_samples).raw_dir};
wt_indices = find(strcmp(genotypes, 'wt'));
bynMyo1c_indices = find(strcmp(genotypes, 'bynMyo1c'));
refIDs = {2, 6};

%% 1. radius kymographs, interpolated between folds
fprintf('\n========== PHASE 1: Radius-based timestamping ==========\n');
Radii_interp = cell(numel(dset_samples), 1);
for dsetID = 1:numel(dset_samples)
    s = dset_samples(dsetID);
    tubi = load_tubi(P, s);
    n = P.samples(s).n_t_timealign;
    if ~isempty(n), tubi.xp.fileMeta.timePoints = tubi.xp.fileMeta.timePoints(1:n); end

    % folds straight from TubULAR's fold file (tubi.getFeatures('folds') reads
    % the wrong field name in this TubULAR version)
    foldfn = tubi.fileName.features.folds;
    if ~exist(foldfn, 'file')
        error('Run tubi.identifyFolds first for sample %d (missing %s)', s, foldfn);
    end
    load(foldfn, 'folds');
    fold1_norm = double(folds(:, 1)) / tubi.nU;
    fold2_norm = double(folds(:, 2)) / tubi.nU;
    fold3_norm = double(folds(:, 3)) / tubi.nU;

    radii = zeros(length(tubi.xp.fileMeta.timePoints), tubi.nU);
    for tidx = 1:length(tubi.xp.fileMeta.timePoints)
        tubi.setTime(tubi.xp.fileMeta.timePoints(tidx))
        mesh = tubi.getCurrentSPCutMesh();
        radii(tidx,:) = mean(mesh.radii_from_mean_uniform_rs, 2);
    end
    Radii_interp{dsetID} = interpolate_spatially_between_folds(radii, fold1_norm, fold2_norm, fold3_norm, ...
        'soverL', linspace(0, 1, tubi.nU), 'chamber_sizes', [100, 99, 99, 99]);
    fprintf('sample %d (dataset %d): radii %s\n', s, dsetID, mat2str(size(Radii_interp{dsetID})));
end

%% 2. shifts within each genotype
rel_shift_fn = P.timeshifts;
if exist(rel_shift_fn, 'file') && ~overwrite
    error('%s exists; set overwrite = true to recompute', rel_shift_fn);
end
[timeshifts, optimal_shifts, correlations, shift_values] = ...
    compute_radius_timestamps(Radii_interp, wt_indices, bynMyo1c_indices, outdir, refIDs);
relative_shifts = struct();
genotype_names = {'wt', 'bynMyo1c'};
for g = 1:length(timeshifts)
    if ~isempty(timeshifts{g})
        relative_shifts.(genotype_names{g}) = timeshifts{g};
    end
end

readme = struct();
readme.description = 'Radius-based time shift alignment between datasets';
readme.relative_shifts = 'Struct with fields wt and bynMyo1c containing time shifts relative to reference dataset for each genotype';
readme.timeshifts = 'Cell array {genotype}(dataset_position) of optimal time shifts (in timepoints) for each dataset relative to its genotype reference';
readme.optimal_shifts = 'Cell array {genotype}{dataset_position} of computed optimal shifts from cross-correlation';
readme.dataDirs = 'Cell array of data directory paths for all datasets';
readme.genotypes = 'Cell array of genotype labels (''wt'' or ''bynMyo1c'') for each dataset';
readme.wt_indices = 'Indices of WT datasets in dataDirs';
readme.bynMyo1c_indices = 'Indices of bynMyo1c datasets in dataDirs';
readme.refIDs = 'Cell array {genotype_idx} of reference dataset IDs used for alignment within each genotype';
readme.cross_genotype_shift = 'Optimal time shift (in timepoints) of bynMyo1c reference relative to WT reference';
readme.wt_ref_idx = 'Dataset index of WT reference used for cross-genotype comparison';
readme.byn_ref_idx = 'Dataset index of bynMyo1c reference used for cross-genotype comparison';
readme.max_corr_cross = 'Maximum correlation value at optimal cross-genotype shift';
readme.corr_vals_cross = 'Array of correlation values at all tested shifts for cross-genotype comparison';
readme.shift_vals_cross = 'Array of time shift values tested for cross-genotype comparison';
save(rel_shift_fn, 'relative_shifts', 'timeshifts', 'optimal_shifts', ...
    'dataDirs', 'genotypes', 'wt_indices', 'bynMyo1c_indices', 'refIDs', 'readme');

%% 3. cross-genotype shift between the two references
wt_ref_idx = refIDs{1};
byn_ref_idx = refIDs{2};
[cross_genotype_shift, max_corr_cross, corr_vals_cross, shift_vals_cross] = ...
    compute_cross_genotype_timing(Radii_interp, wt_ref_idx, byn_ref_idx, wt_indices, bynMyo1c_indices, outdir);
save(rel_shift_fn, 'cross_genotype_shift', 'wt_ref_idx', 'byn_ref_idx', ...
    'max_corr_cross', 'corr_vals_cross', 'shift_vals_cross', 'readme', '-append');
fprintf('wrote %s\n', rel_shift_fn);
