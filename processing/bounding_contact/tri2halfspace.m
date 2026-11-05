function [n, d] = tri2halfspace(V1, V2, V3, centroid)
% returns normal vector and float parameter (n, d) such that if dot(n, V) < d, V is inside the halfspace.
% uses centroid to fix the sign of the inequality.

n = cross(V2-V1, V3-V2, 2);
n = n ./ vecnorm(n,2,2);
d = dot(n, V1, 2);
for i = 1:size(n,1) % n is Nx3 array
    if dot(centroid, n(i,:)) > d(i)
        n(i,:) = -n(i,:);
        d(i) = -d(i);
    end
end


end