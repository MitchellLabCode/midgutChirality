function hs = hoop_arclength(data)
% hoop_arclength  Arclength fraction of each TubULAR mesh hoop, per frame.
%
%   hs = hoop_arclength(data)
%
% Returns hs (n_frames x n_hoops): cumulative arclength of the hoop
% centroids, normalised to [0, 1]. Use it to put per-hoop quantities (e.g.
% data.radii) on an arclength axis, and to place folds given as hoop
% indices: hs(t, data.folds(t,k)). Hoops are not evenly spaced in arclength,
% so hoop_index/n_hoops is not s/L.

nt = numel(data.meshes);
v1 = data.meshes{1}.v;
n_hoop = 100;                       % nU of the mesh grid
n_v = size(v1, 1) / n_hoop;
hs = zeros(nt, n_hoop);
for t = 1:nt
    V  = reshape(data.meshes{t}.v, [n_hoop, n_v, 3]);
    hc = squeeze(mean(V, 2));                                % hoop centroids
    ss = [0; cumsum(vecnorm(diff(hc, 1, 1), 2, 2))];
    hs(t, :) = ss / ss(end);
end
end
