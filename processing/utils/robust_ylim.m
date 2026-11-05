function yl = robust_ylim(xl, varargin)
% robust_ylim  Shared [0, y_top] limits for the time-aligned summary plots.
%
%   yl = robust_ylim(xl, t1, y1, t2, y2, ...)
%
% Pools every single-embryo trace value y_k whose aligned time t_k falls
% inside the plotted window xl = [t_min, t_max], and returns
% [0, 1.05 * 99th percentile]. The ensemble means and the bulk of the
% individual traces stay visible, while the sparsely-sampled tails of the
% aligned time axis (one or two embryos per bin) do not set the scale.
% Traces are pooled across all panels so the limits can be shared.
%
% Also prints, per series, the fraction of in-window points above y_top, so
% it is visible which panels are being clipped.

pool = [];
ys = cell(1, numel(varargin)/2);
for k = 1:2:numel(varargin)
    t = varargin{k}(:);
    y = varargin{k+1}(:);
    y = y(t >= xl(1) & t <= xl(2));
    ys{(k+1)/2} = y(isfinite(y));
    pool = [pool; y]; %#ok<AGROW>
end
pool = sort(pool(isfinite(pool)));
yl = [0, 1.05 * pool(ceil(0.99 * numel(pool)))];

% clipping report
msg = sprintf('robust_ylim: y_top = %.3g, clipped:', yl(2));
for k = 1:numel(ys)
    name = inputname(2*k + 1);   % y-argument of pair k (arg 1 is xl)
    if isempty(name); name = sprintf('series%d', k); end
    msg = [msg, sprintf(' %s %.1f%%', name, 100*mean(ys{k} > yl(2)))]; %#ok<AGROW>
end
fprintf('%s\n', msg);
end
