function [tri, X, Y, Z] = spherical_patch(bp)
% spherical_patch  Triangulated region of the unit sphere enclosed by a closed
% loop of unit vectors, for drawing the signed area q (qmeasure.m).
%
%   [tri, X, Y, Z] = spherical_patch(bp)
%   trisurf(tri, 1.01*X, 1.01*Y, 1.01*Z, ...)
%
% bp : (n x 3) closed loop of unit vectors (first point repeated at the end)
%
% The loop is mapped to spherical angles (theta, phi), a regular angle grid
% is clipped to the enclosed polygon and triangulated with an alpha shape.
% Assumes the loop does not cross phi = +-pi or enclose a pole.

theta_b = acos(bp(:,3));
phi_b   = atan2(bp(:,2), bp(:,1));
[TH, PH] = meshgrid(linspace(0, pi, 400), linspace(-pi, pi, 799));
in = inpolygon(TH(:), PH(:), theta_b, phi_b);
th = TH(in); ph = PH(in);
tri = alphaTriangulation(alphaShape(ph, th));
X = cos(ph) .* sin(th);
Y = sin(ph) .* sin(th);
Z = cos(th);
end
