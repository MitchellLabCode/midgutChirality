function Rc = contact_radius(d)
% contact_radius  Contact radius at every centerline point and time [um]:
% 1.2 x the distance from the centerline point to the nearest surface-mesh
% vertex, i.e. the radius of the sphere centred on the centerline that
% touches the surface, inflated by 20%. Used as the confinement of the rod
% simulations (cosserat/).
%
%   Rc = contact_radius(d)      % d: one element of all_data; Rc (n_t x N)

x = d.cent;
[num_t, num_s, ~] = size(x);
Rc = zeros(num_t, num_s);
for t = 1:min(num_t, numel(d.meshes))
    v = d.meshes{t}.v;
    for s = 1:num_s
        Rc(t,s) = 1.2 * min(vecnorm(v - squeeze(x(t,s,:))', 2, 2));
    end
end
end
