function data_interp = interpolate_spatially_between_folds(data, fold1_norm, fold2_norm, fold3_norm, varargin)
    % Interpolate spatial data (e.g., radii, strains) between fold boundaries
    % Similar to the interpolation in segment_strains_updated.m
    %
    % Usage:
    %   data_interp = interpolate_spatially_between_folds(data, fold1_norm, fold2_norm, fold3_norm)
    %   data_interp = interpolate_spatially_between_folds(data, fold1_norm, fold2_norm, fold3_norm, ...
    %       'soverL', soverL, 'chamber_sizes', [100, 99, 99, 99])
    %
    % Input:
    %   data: matrix of size (n_timepoints x n_spatial_points)
    %   fold1_norm, fold2_norm, fold3_norm: vectors of length n_timepoints with normalized fold positions [0,1]
    %
    % Optional name-value pairs:
    %   'soverL': spatial coordinates of original data points (default: linspace(0, 1, n_spatial_points))
    %   'chamber_sizes': array of 4 integers for number of points in each chamber (default: [100, 99, 99, 99])
    %
    % Output:
    %   data_interp: interpolated matrix of size (n_timepoints x n_interp_pts)
    %               where n_interp_pts = sum(chamber_sizes)
    
    % Parse inputs
    p = inputParser;
    n_spatial = size(data, 2);
    addParameter(p, 'soverL', linspace(0, 1, n_spatial), @isvector);
    addParameter(p, 'chamber_sizes', [100, 99, 99, 99], @isvector);
    parse(p, varargin{:});
    
    soverL = p.Results.soverL(:);  % Ensure column vector
    chamber_sizes = p.Results.chamber_sizes(:);  % Ensure column vector
    
    % Validate inputs
    if size(data, 1) ~= length(fold1_norm) || ...
       length(fold1_norm) ~= length(fold2_norm) || ...
       length(fold2_norm) ~= length(fold3_norm)
        error('Fold vectors must have same length as number of timepoints');
    end
    
    if length(soverL) ~= size(data, 2)
        error('soverL length must match number of spatial points in data');
    end
    
    if length(chamber_sizes) ~= 4
        error('chamber_sizes must have exactly 4 elements (one per chamber)');
    end
    
    % Initialize output
    n_timepoints = size(data, 1);
    n_interp_pts = sum(chamber_sizes);
    data_interp = zeros(n_timepoints, n_interp_pts);
    
    % Process each timepoint
    for t = 1:n_timepoints
        % Get fold positions for this timepoint
        f1 = fold1_norm(t);
        f2 = fold2_norm(t);
        f3 = fold3_norm(t);
        
        % Define chamber boundaries: [0, f1], (f1, f2], (f2, f3], (f3, 1]
        chamber_bounds = [0, f1, f2, f3, 1];
        
        % Process each chamber
        col_idx = 1;
        for chamber_idx = 1:4
            x_start = chamber_bounds(chamber_idx);
            x_end = chamber_bounds(chamber_idx + 1);
            n_pts = chamber_sizes(chamber_idx);
            
            % Find all original points in this chamber
            mask = soverL >= x_start & soverL <= x_end;
            x_orig = soverL(mask);
            
            if sum(mask) > 0
                % Create bin edges for the new grid in this chamber
                bin_edges = linspace(x_start, x_end, n_pts + 1);
                
                % Bin and average the data
                data_binned = zeros(1, n_pts);
                for bin_idx = 1:n_pts
                    bin_mask = x_orig >= bin_edges(bin_idx) & x_orig <= bin_edges(bin_idx + 1);
                    if sum(bin_mask) > 0
                        orig_indices = find(mask);
                        bin_orig_indices = orig_indices(bin_mask);
                        data_binned(bin_idx) = mean(data(t, bin_orig_indices));
                    else
                        % Interpolate if bin is empty, using extrap to handle edge cases
                        bin_center = bin_edges(bin_idx) + (bin_edges(bin_idx+1) - bin_edges(bin_idx))/2;
                        data_binned(bin_idx) = interp1(x_orig, data(t, mask), bin_center, 'linear', 'extrap');
                    end
                end
                
                % Store in output array
                data_interp(t, col_idx:col_idx+n_pts-1) = data_binned;
                col_idx = col_idx + n_pts;
            end
        end
    end
end
