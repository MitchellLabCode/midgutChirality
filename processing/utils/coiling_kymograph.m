function [Q, s_grid] = coiling_kymograph(curves, w, n_skip_start, n_skip_end)
% coiling_kymograph  Quasi-local coiling (qmeasure.m) measure along a centerline,
% for every frame, on a common arclength axis.
%
%   [Q, s_grid] = coiling_kymograph(curves, w, n_skip_start, n_skip_end)
%
% curves       : n_t x N x 3 array, or a 1 x n_t cell of N_t x 3 curves
% w            : half-width of the sliding window, in tangents (window 2w+1)
% n_skip_start : points dropped at the anterior end before computing tangents
% n_skip_end   : points dropped at the posterior end
%
% Q      : n_s x n_t measure, columns = frames
% s_grid : 1 x n_s positions as arclength fractions of the FULL curve, so
%          the measure lines up with fold positions from fold_arclength /
%          hoop arclength, and with s/L of any other kymograph.
%
% Each value sits at the arclength midpoint of the tangent its window is
% centred on. Frames are interpolated onto one grid spanning the range all
% frames cover. Earlier inline versions placed the measure at the index it
% had after dropping n_skip_start points (a ~9% anterior shift for
% n_skip_start = 9) and zero-padded the posterior end.

if iscell(curves)
    n_t = numel(curves);
    get_curve = @(t) curves{t};
else
    n_t = size(curves, 1);
    get_curve = @(t) squeeze(curves(t, :, :));
end

q_t = cell(1, n_t);
s_t = cell(1, n_t);
for t = 1:n_t
    [q_t{t}, s_t{t}] = coiling_measure(get_curve(t), w, n_skip_start, n_skip_end);
end

s_lo = max(cellfun(@(s) s(1),   s_t));
s_hi = min(cellfun(@(s) s(end), s_t));
n_s  = min(cellfun(@numel, q_t));
s_grid = linspace(s_lo, s_hi, n_s);

Q = zeros(n_s, n_t);
for t = 1:n_t
    Q(:, t) = interp1(s_t{t}, q_t{t}, s_grid, 'linear');
end
end

function [q, s_c] = coiling_measure(curve, w, n_skip_start, n_skip_end)
% measure for one frame, and window-centre positions as fractions of the
% full curve's arclength
ss  = [0; cumsum(vecnorm(diff(curve, 1, 1), 2, 2))];
sub = curve(n_skip_start+1 : end-n_skip_end, :);
tangents = diff(sub, [], 1);
tangents = tangents ./ vecnorm(tangents, 2, 2);
n_tan = size(tangents, 1);

centres = (1+w):(n_tan-w);            % tangent index each window is centred on
q   = zeros(numel(centres), 1);
s_c = zeros(1, numel(centres));
for j = 1:numel(centres)
    idx = centres(j);
    [q(j), ~] = qmeasure(tangents(idx-w:idx+w, :));
    p = n_skip_start + idx;           % tangent idx joins full-curve points p, p+1
    s_c(j) = 0.5 * (ss(p) + ss(p+1)) / ss(end);
end
end
