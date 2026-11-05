function [n,d] = boundingplanes(x, y, z)
% finds the convex hull of the list of N points described by x,y,z
% coordinates (all Nx1 arrays)  and returns the normals and parameters of
% the half-space representation of the convex hull polyhedron.
%
%

[k2,av2] = convhull(x, y, z,'Simplify',true);

verts = [x, y, z];

centroid = mean(verts, 1);

V1 = verts(k2(:, 1),:);

V2 = verts(k2(:, 2),:);

V3 = verts(k2(:, 3),:);

[n,d] = tri2halfspace(V1, V2, V3, centroid);


end