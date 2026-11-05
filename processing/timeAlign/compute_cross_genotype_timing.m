%% COMPUTE CROSS-GENOTYPE TIMING
% Performs cross-correlation analysis between reference datasets of different
% genotypes to find optimal time shift.
%
% Inputs:
%   Radii_interp     - Cell array of interpolated radius kymographs (one per dataset)
%   wt_ref_idx       - Index of WT reference dataset
%   byn_ref_idx      - Index of bynMyo1c reference dataset
%   wt_indices       - Indices of all WT datasets
%   bynMyo1c_indices - Indices of all bynMyo1c datasets
%   outdir           - Output directory for saving results
%
% Outputs:
%   cross_genotype_shift - Optimal time shift (bynMyo1c relative to WT)
%   max_corr_cross       - Maximum correlation at optimal shift
%   corr_vals_cross      - Correlation values at all tested shifts
%   shift_vals_cross     - Shift values tested

function [cross_genotype_shift, max_corr_cross, corr_vals_cross, shift_vals_cross] = ...
    compute_cross_genotype_timing(Radii_interp, wt_ref_idx, byn_ref_idx, wt_indices, bynMyo1c_indices, outdir)

    fprintf('\n========== Cross-genotype timing: WT ref vs bynMyo1c ref ==========\n');
    
    refRadii_wt = Radii_interp{wt_ref_idx};
    refRadii_byn = Radii_interp{byn_ref_idx};

    fprintf('Comparing WT reference (dataset %d) with bynMyo1c reference (dataset %d)\n', ...
        wt_ref_idx, byn_ref_idx);

    % Compute correlations with time shifts
    options = struct();
    options.debug = false;
    options.max_shift = 30;
    [corr_vals_cross, shift_vals_cross, optimal_shift_cross, ~, ~, max_corr_cross, corrsm_cross] = ...
        compute_kymoshifted_temporal_xcorrelation(refRadii_byn, refRadii_wt, options);

    fprintf('Optimal time shift (bynMyo1c ref relative to WT ref): %d timepoints (correlation = %.4f)\n', ...
        optimal_shift_cross, max_corr_cross);

    % Visualize cross-genotype comparison
    n_tp_wt = size(refRadii_wt, 1);
    n_tp_byn = size(refRadii_byn, 1);
    max_tp_cross = max(n_tp_wt, n_tp_byn + abs(optimal_shift_cross));
    n_spatial = size(refRadii_byn, 2);

    % Create shifted bynMyo1c radii
    refRadii_byn_shifted = NaN(max_tp_cross, n_spatial);
    if optimal_shift_cross >= 0
        refRadii_byn_shifted(optimal_shift_cross+1:optimal_shift_cross+n_tp_byn, :) = refRadii_byn;
    else
        refRadii_byn_shifted(1:n_tp_byn+optimal_shift_cross, :) = refRadii_byn(-optimal_shift_cross+1:end, :);
    end

    close all
    figure('Name', 'Cross-Genotype Radii Comparison', 'Position', [100, 100, 1600, 400]);

    subplot(1, 4, 1)
    imagesc(refRadii_wt)
    title(sprintf('WT Reference Radii\nWT #%d', find(wt_indices == wt_ref_idx)))
    colorbar
    set(gca, 'YLim', [0.5, max_tp_cross + 0.5]);

    subplot(1, 4, 2)
    imagesc(refRadii_byn)
    title(sprintf('bynMyo1c Reference Radii (unshifted)\nbyn>Myo1C #%d', find(bynMyo1c_indices == byn_ref_idx)))
    colorbar
    set(gca, 'YLim', [0.5, max_tp_cross + 0.5]);

    subplot(1, 4, 3)
    imagesc(refRadii_byn_shifted)
    title(sprintf('bynMyo1c Reference Radii (shifted by %d tp)\nbyn>Myo1C #%d', optimal_shift_cross, find(bynMyo1c_indices == byn_ref_idx)))
    colorbar
    set(gca, 'YLim', [0.5, max_tp_cross + 0.5]);

    subplot(1, 4, 4)
    plot(shift_vals_cross, corr_vals_cross, 'o-', 'LineWidth', 2);
    hold on;
    plot(shift_vals_cross, corrsm_cross, '-', 'LineWidth', 1);
    plot(optimal_shift_cross * [1,1], max_corr_cross * [0,1], 'Color', 'r', 'LineStyle', '--');
    xlabel('Time Shift (timepoints)');
    ylabel('Mean Correlation');
    title(sprintf('Cross-Genotype Cross-Correlation\nWT #%d vs byn>Myo1C #%d', ...
        find(wt_indices == wt_ref_idx), find(bynMyo1c_indices == byn_ref_idx)));
    grid on;

    % Return cross-genotype timing result
    cross_genotype_shift = optimal_shift_cross;

    saveas(gcf, fullfile(outdir, 'cross_genotype_Radii_CrossCorrelation.png'));
    if usejava('desktop')            % interactive session: inspect before continuing
        disp('close fig to continue')
        waitfor(gcf);
    end

end
