function [t_mid, v_mean] = bin_means(t, v, n_bins)
% bin_means  Mean of pooled (time, value) samples in n_bins equal time bins
% spanning the pooled range; empty bins are dropped.
%
%   [t_mid, v_mean] = bin_means(t, v, n_bins)

t = t(:); v = v(:);
edges = linspace(min(t), max(t), n_bins + 1);
[~, ~, loc] = histcounts(t, edges);
valid = loc > 0;
v_mean = accumarray(loc(valid), v(valid)) ./ accumarray(loc(valid), 1);
t_mid = 0.5 * (edges(1:end-1) + edges(2:end));
has_data = false(n_bins, 1);
has_data(unique(loc(valid))) = true;
t_mid  = t_mid(has_data);
v_mean = v_mean(has_data);
end
