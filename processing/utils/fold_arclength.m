function fs = fold_arclength(data, coord)
% fold_arclength  Fold positions as arclength fractions, for kymograph overlays.
%
%   fs = fold_arclength(data, coord)
%
% coord = 'x'    : fraction along the full centerline data.x, i.e.
%                  TubULAR's own fold arclength ssfold ./ ssmax
%         'cent' : fraction along data.cent, which drops the first
%                  n_cut = size(data.x,2) - size(data.cent,2) points of data.x
%
% Returns fs with the size of data.ssfold (one row per frame, one column per
% fold), NaN where undefined. Rows line up with data.folds, so the usual
% onset masks (folds(:,k) > 0) still select the valid rows.
%
% Replaces double(folds)/num_s, which treats the TubULAR hoop index as if
% hoops were evenly spaced in arclength. They are not: index/num_s put the
% folds up to ~0.1 of the gut length too far posterior.

fs = data.ssfold ./ data.ssmax;
fs(~isfinite(fs) | fs < 0 | fs > 1) = NaN;

switch coord
    case 'x'
    case 'cent'
        n_cut = size(data.x, 2) - size(data.cent, 2);
        nt = min(size(fs, 1), size(data.x, 1));
        for t = 1:nt
            xt = squeeze(data.x(t, :, :));
            ss = [0; cumsum(vecnorm(diff(xt), 2, 2))];
            c  = ss(n_cut + 1) / ss(end);      % arclength fraction cut off
            fs(t, :) = (fs(t, :) - c) / (1 - c);
        end
    otherwise
        error('fold_arclength: coord must be ''x'' or ''cent''');
end
end
