function [q, centres] = coiling_profile(curv, w)
% coiling_profile  Quasi-local chirality q along one curve (qmeasure.m over a
% window of 2w+1 unit tangents), for the 3D renders and movies.
%
%   [q, centres] = coiling_profile(curv, w)
%
% curv    : (N x 3) curve
% q       : value at each window centre (unsmoothed)
% centres : index of the tangent (and curve point) each window is centred on;
%           plot q at curv(centres, :)

tangents = diff(curv, [], 1);
tangents = tangents ./ vecnorm(tangents, 2, 2);
centres = (1+w):(size(curv, 1) - w - 1);
q = zeros(numel(centres), 1);
for m = 1:numel(centres)
    j = centres(m);
    q(m) = qmeasure(tangents(j-w:j+w, :));
end
end
