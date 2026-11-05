%% COMPUTE RADIUS-BASED TIMESTAMPS
% Performs cross-correlation analysis on interpolated radii to find optimal
% time shifts between datasets within each genotype.
%
% Inputs:
%   Radii_interp   - Cell array of interpolated radius kymographs (one per dataset)
%   wt_indices     - Indices of WT datasets
%   bynMyo1c_indices - Indices of bynMyo1c datasets
%   outdir         - Output directory for saving results
%   refIDs         - [optional] Reference dataset IDs for each genotype (default: {3, 4})
%
% Outputs:
%   timeshifts     - Cell array of optimal time shifts for each genotype
%   optimal_shifts - Cell array of all optimal shifts computed
%   correlations   - Cell array of correlation curves
%   shift_values   - Cell array of shift values tested
%
% Called from timelag_fromRadii.m after radii interpolation.

function [timeshifts, optimal_shifts, correlations, shift_values] = ...
    compute_radius_timestamps(Radii_interp, wt_indices, bynMyo1c_indices, outdir, refIDs)

    if nargin < 5
        refIDs = {3, 4}; % Default reference dataset IDs (one per genotype)
    end

    % Define genotypes and their corresponding indices
    genotype_list = {'wt', 'bynMyo1c'};
    indices_list = {wt_indices, bynMyo1c_indices};

    % Store results
    optimal_shifts = {};
    correlations = {};
    shift_values = {};
    timeshifts = cell(length(genotype_list), 1);

    for g = 1:length(genotype_list)
        genotype = genotype_list{g};
        indices = indices_list{g};
        
        if length(indices) >= 2
            fprintf('\n=== Cross-correlation analysis for %s radii ===\n', genotype);
            
            % Find the position of the reference dataset in indices
            ref_pos = find(indices == refIDs{g});
            if isempty(ref_pos)
                error('Reference dataset %d not found in %s genotype indices %s', ...
                    refIDs{g}, genotype, mat2str(indices));
            end
            
            % Use the reference dataset
            dset1_idx = indices(ref_pos);
            refRadii = Radii_interp{dset1_idx};
            
            % Compare all other positions to the reference
            datasets2do = setdiff(1:length(indices), ref_pos);
            fprintf('Reference dataset: %d\n', dset1_idx);
            fprintf('Comparing datasets: %s\n', mat2str(indices(datasets2do)));
            disp('datasets to do: ')
            disp(datasets2do)
            
            for compareIdx = datasets2do
                dset2_idx = indices(compareIdx);
                radii2 = Radii_interp{dset2_idx};
                
                fprintf('Comparing %s dataset %d (%d timepoints) with dataset %d (%d timepoints)\n', ...
                    genotype, dset1_idx, size(refRadii, 1), dset2_idx, size(radii2, 1));
                
                % Compute correlations with time shifts
                options = struct();
                options.debug = false;
                options.max_shift = 30;
                [corr_vals, shift_vals, ...
                    optimal_shift, final_overlap1, final_overlap2, max_corr, corrsm] =...
                    compute_kymoshifted_temporal_xcorrelation(radii2, refRadii, options);
                
                % Store results
                optimal_shifts{g}{compareIdx} = optimal_shift;
                correlations{g}{compareIdx} = corr_vals;
                shift_values{g}{compareIdx} = shift_vals;
                
                fprintf('Optimal time shift: %d timepoints (correlation = %.4f)\n', optimal_shift, max_corr);
                fprintf('All shifts and correlations:\n');
                disp(table(shift_vals(:), corr_vals(:), 'VariableNames', {'TimeShift', 'Correlation'}));
                
                % Plot
                close all
                figure('Name', sprintf('%s Radii Cross-Correlation vs Time Shift', genotype), 'Position', [100, 100, 1600, 400]);
                
                % Get consistent time axis limits for aligned y-axis
                n_tp_ref = size(refRadii, 1);
                n_tp_comp = size(radii2, 1);
                % Extended size to accommodate shift without truncation
                max_tp = max(n_tp_ref, n_tp_comp + abs(optimal_shift));
                
                % Determine genotype labels with indices
                if strcmp(genotype, 'wt')
                    ref_genotype_label = sprintf('WT #%d (ref)', find(wt_indices == dset1_idx));
                    comp_genotype_label = sprintf('WT #%d', find(wt_indices == dset2_idx));
                else
                    ref_genotype_label = sprintf('byn>Myo1C #%d (ref)', find(bynMyo1c_indices == dset1_idx));
                    comp_genotype_label = sprintf('byn>Myo1C #%d', find(bynMyo1c_indices == dset2_idx));
                end
                
                subplot(1, 4, 1)
                imagesc(refRadii)
                title(sprintf('Reference Radii\n%s', ref_genotype_label))
                colorbar
                set(gca, 'YLim', [0.5, max_tp + 0.5]);
                
                subplot(1, 4, 2)
                imagesc(radii2)
                title(sprintf('Comparison Radii (unshifted)\n%s', comp_genotype_label))
                colorbar
                set(gca, 'YLim', [0.5, max_tp + 0.5]);
                
                % Create shifted radii for visualization (linear shift, no wrapping)
                % Extend array size to accommodate shift without truncation
                n_spatial = size(radii2, 2);
                radii2_shifted = NaN(max_tp, n_spatial);
                if optimal_shift >= 0
                    % Positive shift: move data down, pad top with NaN
                    radii2_shifted(optimal_shift+1:optimal_shift+n_tp_comp, :) = radii2;
                else
                    % Negative shift: move data up, pad bottom with NaN
                    shift_start = 1 - optimal_shift;
                    radii2_shifted(1:n_tp_comp+optimal_shift, :) = radii2(-optimal_shift+1:end, :);
                end
                subplot(1, 4, 3)
                imagesc(radii2_shifted)
                title(sprintf('Comparison Radii (shifted by %d tp)\n%s', optimal_shift, comp_genotype_label))
                colorbar
                set(gca, 'YLim', [0.5, max_tp + 0.5]);
                
                subplot(1, 4, 4)
                plot(shift_vals, corr_vals, 'o-', 'LineWidth', 2);
                hold on;
                plot(shift_vals, corrsm, '-', 'LineWidth', 1) ;
                plot(optimal_shift * [1,1], max_corr * [0,1], 'Color', 'r', 'LineStyle', '--');
                
                xlabel('Time Shift (timepoints)');
                ylabel('Mean Correlation');
                title(sprintf('%s Cross-Correlation\n%s vs %s', genotype, ref_genotype_label, comp_genotype_label));
                grid on;

                set(gcf, 'color', 'w');
                % Save figure
                fig_fn = fullfile(outdir, sprintf('%s_Radii_CrossCorrelation_Dataset%d_vs_Dataset%d.png', ...
                    genotype, dset1_idx, dset2_idx));
                saveas(gcf, fig_fn);

                % disp('close fig to continue')
                % waitfor(gcf);
        
                timeshifts{g}{compareIdx} = optimal_shift;
            end

        else
            fprintf('Skipping %s: only %d dataset(s) available (need at least 2)\n', genotype, length(indices));
        end

        % Save the timestamping/timeshifts
        ofn = fullfile(outdir, 'timeshifts_radii.mat');
        save(ofn, 'optimal_shifts', 'timeshifts', 'shift_values', 'correlations')
    end

end
