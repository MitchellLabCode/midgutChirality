function [sum_area, areas] = qmeasure(vectors)
% QMEASURE  Quasi-local chirality measure q: signed area on the unit sphere
% swept by a sequence of unit tangent vectors (spherical-triangle sum;
% formerly stare.m, Spherical Triangle Area Rotation Estimator).
% input:
% - vectors Nx3: a set of unit vectors
% Output:
% - the signed area of the (perhaps non-convex) spherical polygon defined
% by the unit vectors

function A=SArea(a, b, c)
    A = 2 * atan( dot(cross(a,b), c) / (1 + dot(a,b) + dot(b,c) + dot(a,c)) );
end

num_vectors = size(vectors, 1);
if num_vectors <3
    areas = 0;
    sum_area = 0;
    disp('Warning: not enough points. returned 0 area')
else
    areas = zeros(num_vectors-2,1);
    for k = 2:(num_vectors-1)
        areas(k-1) = SArea(vectors(1,:), vectors(k,:), vectors(k+1,:));
    end
    sum_area = sum(areas);
end % else
end %function
