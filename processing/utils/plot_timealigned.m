function fig = plot_timealigned(S, rows, opts)
% plot_timealigned  Time-aligned summary figure: one row per set of measures,
% one column per region (anterior, posterior, whole gut); thin lines are the
% individual embryos (WT solid, byn>Myo1C dashed), thick lines the binned
% mean over embryos (bin_means.m). Each row has its own y-limits, shared by
% its three panels (robust_ylim.m).
%
%   fig = plot_timealigned(S, rows, opts)
%
% S    : 1 x n_samples struct array with S(i).is_WT and one measure group per
%        field (e.g. S(i).E = compute_deformation_measures(...).E); a group has
%        a time vector t (n x 1) and (n x 3) measures, columns [A P all]
% rows : struct array, one per figure row, with fields
%          src     name of the measure group in S ('E', 'omega', 'mu')
%          names   cell of measure names in that group
%          colors  cell of RGB triplets, one per name
%          labels  cell of legend labels (legend in the row's first panel)
%          ylabel  y-axis label
%          ylim_from  'all' (traces of all three panels) or 'whole'
%                     (whole-gut traces only, e.g. for integrals)
% opts : optional, fields n_bins (30), xl ([-50 150]), position

if nargin < 3, opts = struct(); end
if ~isfield(opts, 'n_bins'), opts.n_bins = 30; end
if ~isfield(opts, 'xl'), opts.xl = [-50, 150]; end
nr = numel(rows);
heights = [420, 820, 1200];
if ~isfield(opts, 'position'), opts.position = [100, 100, 1250, heights(min(nr, 3))]; end
regions = {'anterior', 'posterior', 'whole gut'};

fig = figure('Color', 'w', 'Position', opts.position);
for r = 1:nr
    R = rows(r);
    if ~isfield(R, 'ylim_from') || isempty(R.ylim_from), R.ylim_from = 'all'; end
    G = arrayfun(@(s) s.(R.src), S);                       % this row's measure group
    T = vertcat(G.t);
    pairs = {};
    for c = 1:3
        ax = subplot(nr, 3, 3*(r-1) + c); hold(ax, 'on')
        for i = 1:numel(S)                                 % individual embryos
            ls = '-'; if ~S(i).is_WT, ls = '--'; end
            for m = 1:numel(R.names)
                plot(ax, G(i).t, G(i).(R.names{m})(:,c), 'Color', R.colors{m}, ...
                     'LineStyle', ls, 'LineWidth', 1, 'HandleVisibility', 'off')
            end
        end
        h = gobjects(1, numel(R.names));
        for m = 1:numel(R.names)                           % binned mean over embryos
            V = cell2mat(arrayfun(@(g) g.(R.names{m})(:,c), G(:), 'UniformOutput', false));
            [tm, vm] = bin_means(T, V, opts.n_bins);
            h(m) = plot(ax, tm, vm, 'Color', R.colors{m}, 'LineWidth', 4);
            if strcmp(R.ylim_from, 'all') || c == 3
                pairs = [pairs, {T, V}]; %#ok<AGROW>
            end
        end
        if r == 1, title(ax, regions{c}); end
        if c == 1
            ylabel(ax, R.ylabel)
            if r == 1 || ~isequal(R.labels, rows(r-1).labels)   % legend where the labels change
                legend(h, R.labels, 'Location', 'northwest'); legend(ax, 'boxoff')
            end
        end
        if r == nr, xlabel(ax, 'Time [min, aligned]'); end
        box(ax, 'on'); axis(ax, 'square')
    end
    yl = robust_ylim(opts.xl, pairs{:});
    for c = 1:3
        ax = subplot(nr, 3, 3*(r-1) + c);
        ylim(ax, yl); xlim(ax, opts.xl)
    end
end
end
