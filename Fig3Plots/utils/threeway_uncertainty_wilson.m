function out = threeway_uncertainty_wilson(x, y, z)
% THREEWAY_UNCERTAINTY_WILSON
%   Wilson score 95% confidence intervals for category proportions
%   in a single 3-outcome experiment with counts (x,y,z), n = x+y+z.
%
%   Output (struct):
%     out.n                  : total count
%     out.phat               : [px, py, pz]
%     out.ci95               : 3x2 matrix; rows = {x,y,z}, columns = [L U]
%     out.table              : [x; y; z]
%     out.details.x/y/z      : per-category struct with fields phat, L, U, halfwidth
%
%   Notes (frequentist, one-sample per category):
%     Wilson CI for p = k/n with z = 1.96:
%       den    = 1 + z^2/n
%       center = (phat + z^2/(2n)) / den
%       half   = z * sqrt( phat*(1-phat)/n + z^2/(4n^2) ) / den
%       CI     = [center - half, center + half] clipped to [0,1]
%
%   Example:
%     out = threeway_uncertainty_wilson(23,18,9)

    arguments
        x (1,1) double {mustBeNonnegative, mustBeInteger}
        y (1,1) double {mustBeNonnegative, mustBeInteger}
        z (1,1) double {mustBeNonnegative, mustBeInteger}
    end

    n = x + y + z;
    if n == 0
        error('threeway_uncertainty_wilson:NoData', 'x+y+z must be > 0.');
    end

    % Fixed 95% CI
    alpha = 0.05;
    zcrit = norminv(1 - alpha/2);

    % Helper: Wilson CI for count k out of n
    function S = wilson_one(k, n_, zc)
        ph = k / n_;
        den = 1 + (zc^2)/n_;
        center = (ph + (zc^2)/(2*n_)) / den;
        half   = zc * sqrt( ph*(1-ph)/n_ + (zc^2)/(4*n_^2) ) / den;
        L = max(0, center - half);
        U = min(1, center + half);
        S = struct('phat', ph, 'L', L, 'U', U, 'halfwidth', half);
    end

    Sx = wilson_one(x, n, zcrit);
    Sy = wilson_one(y, n, zcrit);
    Sz = wilson_one(z, n, zcrit);

    out = struct();
    out.n = n;
    out.table = [x; y; z];
    out.phat = [Sx.phat, Sy.phat, Sz.phat];
    out.ci95 = [Sx.L, Sx.U; Sy.L, Sy.U; Sz.L, Sz.U];
    out.details = struct('x', Sx, 'y', Sy, 'z', Sz);
end
