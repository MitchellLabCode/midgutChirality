function [correlations, time_shifts, optimal_shift, final_overlap1, final_overlap2, max_corr, corrsm]... 
        = compute_kymoshifted_temporal_xcorrelation(queryData, refData, options)
    % Compute correlation between two radii heatmaps with time shifts
    % This solves the shift that must be added to the timeline for dataset
    % 1 in order for it to match dataset 2. That is, 
    % queryData(optimal_shift:end) ~= refData(1:ntps), where ntps is the 
    % number of timepoints ranging from optimal_shift to size(queryData, 1)
    %
    % queryData, refData: matrices of size (n_timepoints x n_spatial_pos)
    % options: struct with fields
    %   debug: bool
    %   max_shift: maximum time shift to test (in timepoints)
    % 
    % Returns:
    %   correlations: correlation values for each shift
    %   time_shifts: the shift values tested
    
    if nargin < 3
        options = struct() ;
    end
    if isfield(options, 'max_shift')
        max_shift = options.max_shift ;
    else
        max_shift = size(queryData, 1) / 2;  % Default: half the time range
    end
    
    if isfield(options, 'savgol2')
        savgol2 = options.savgol2 ;
    else
        savgol2 = 7;  % Default: savgol filtering window for computing optimum
    end
    
    if isfield(options, 'debug')
        debug = options.debug ;
    else
        debug = false;  % Default: false
    end
    
    n_t1 = size(queryData, 1);
    n_t2 = size(refData, 1);
    n_s = size(queryData, 2);  % Spatial dimension
    
    % Ensure both have same spatial dimension
    if size(refData, 2) ~= n_s
        error('Radii heatmaps must have same spatial dimension');
    end
    
    % Generate time shifts to test (negative = shift radii2 backward in time)
    time_shifts = -max_shift:max_shift;
    correlations = zeros(size(time_shifts));
    
    % Normalize each row (timepoint) across spatial dimension for correlation
    radii1_norm = queryData - nanmean(queryData, 2);
    radii2_norm = refData - nanmean(refData, 2);
    
    radii1_norm = radii1_norm ./ (std(radii1_norm, [], 2) + eps);
    radii2_norm = radii2_norm ./ (std(radii2_norm, [], 2) + eps);
    
    % For each time shift, compute correlation of overlapping rows
    for shift_idx = 1:length(time_shifts)
        shift = time_shifts(shift_idx);
        
        % Determine overlap region
        if shift >= 0
            % radii2 is shifted forward (later in time), so compare early times of radii1 with late times of radii2
            t1_start = 1;
            t1_end = n_t1 - shift;
            t2_start = 1 + shift;
            t2_end = n_t2;
        else
            % radii2 is shifted backward (earlier in time)
            t1_start = 1 - shift;
            t1_end = n_t1;
            t2_start = 1;
            t2_end = n_t2 + shift;
        end
        
        % Skip if no overlap
        if t1_end < t1_start || t2_end < t2_start
            correlations(shift_idx) = NaN;
            continue;
        end
        
        % Ensure both overlaps have the same size
        overlap_len1 = t1_end - t1_start + 1;
        overlap_len2 = t2_end - t2_start + 1;
        overlap_len = min(overlap_len1, overlap_len2);
        
        if overlap_len <= 0
            correlations(shift_idx) = NaN;
            continue;
        end
        
        % Extract overlapping regions (constrained to same length)
        overlap1 = radii1_norm(t1_start:t1_start+overlap_len-1, :);
        overlap2 = radii2_norm(t2_start:t2_start+overlap_len-1, :);
        
        if debug
            figure(1)
            clf
            subplot(1, 2, 1)
            imagesc(overlap1)
            title('kymo1')
            subplot(1, 2, 2)
            imagesc(overlap2)
            title(['kymo2, shift=' num2str(shift)])
            pause(0.1)
        end

        % Compute mean correlation across all spatial points
        % For each spatial position s, compute correlation between temporal profiles
        corr_per_spatial = zeros(n_s, 1);
        for s = 1:n_s
            corr_per_spatial(s) = corr(overlap1(:, s), overlap2(:, s));
        end
        
        % Mean correlation across spatial dimension
        correlations(shift_idx) = mean(corr_per_spatial, 'omitnan');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REPORTING OPTIMUM
    if nargout > 2
        % Find optimal shift and extract overlap for visualization
        if savgol2 > 0
            corrsm = savgol(correlations, 2, savgol2) ;
        else
            corrsm = correlations ;
        end
        [max_corr, max_idx] = max(corrsm);
        optimal_shift = time_shifts(max_idx);
        
        % Recompute overlap indices for optimal shift
        if optimal_shift >= 0
            t1_start = 1;
            t1_end = n_t1 - optimal_shift;
            t2_start = 1 + optimal_shift;
            t2_end = n_t2;
        else
            t1_start = 1 - optimal_shift;
            t1_end = n_t1;
            t2_start = 1;
            t2_end = n_t2 + optimal_shift;
        end
        
        overlap_len = min(t1_end - t1_start + 1, t2_end - t2_start + 1);
        final_overlap1 = radii1_norm(t1_start:t1_start+overlap_len-1, :);
        final_overlap2 = radii2_norm(t2_start:t2_start+overlap_len-1, :);
        
        fprintf('Optimal shift: %d timepoints (correlation: %.4f)\n', optimal_shift, max_corr);
        
        % Visualize optimal alignment
        if debug
            figure
            subplot(1, 3, 1)
            imagesc(final_overlap1)
            title('Radii 1 overlap')
            xlabel('Spatial position')
            ylabel('Time')
            colorbar
            
            subplot(1, 3, 2)
            imagesc(final_overlap2)
            title(sprintf('Radii 2 overlap (shift=%d)', optimal_shift))
            xlabel('Spatial position')
            ylabel('Time')
            colorbar
    
            subplot(1, 3, 3)
            plot(time_shifts, correlations, '.-')
            hold on;
            plot(time_shifts, corrsm, '-')
            xlabel('Time shifts [frames]')
            ylabel('Correlations')
            hold on;
            plot(optimal_shift * [1,1], max_corr*[0,1], 'Color', 'r', 'LineStyle', '--')
        end
    end
end