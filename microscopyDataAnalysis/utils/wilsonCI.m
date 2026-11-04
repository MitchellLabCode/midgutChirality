function [pHat, lo, hi] = wilsonCI(x, n, alpha)
% WILSONCI  Wilson score confidence interval for a binomial proportion.
%   [pHat, lo, hi] = wilsonCI(x, n, alpha)
%   x     : number of successes
%   n     : number of trials
%   alpha : significance level (default 0.05 for a 95% CI)

    if nargin < 3, alpha = 0.05; end

    % norminv(1 - alpha/2) via the erfinv identity (no Stats Toolbox needed).
    z = sqrt(2) * erfinv(1 - alpha);
    pHat = x / n;

    denom  = 1 + z^2 / n;
    center = (pHat + z^2 / (2*n)) / denom;
    half   = (z / denom) * sqrt(pHat*(1 - pHat)/n + z^2/(4*n^2));

    lo = max(0, center - half);
    hi = min(1, center + half);
end
