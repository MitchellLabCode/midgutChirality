% tubi_to_meshrs.m
% Exports, for every embryo, the TubULAR cut meshes of all timepoints in the
% anatomical (APDV, "RS") frame, and the mean radius of each mesh hoop.
%
% Writes P.meshes/meshes_sample<N>.mat with
%   meshes : (n_t x 1) cell of struct(f, v), v (nU*nV x 3) [um], nU = nV = 100
%   radii  : (n_t x nU) mean distance of each hoop's vertices to its centroid [um]
% Timepoints are cropped as in config_paths.m (n_t_meshes).
% Requires TubULAR. Run before process_all_V6.m.

clear; close all
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'), here); P = config_paths();
if ~exist(P.meshes, 'dir'); mkdir(P.meshes); end

nU = 100; nV = 100;
for idx_sample = 1:numel(P.samples)
    tubi = load_tubi(P, idx_sample);
    timePoints = tubi.xp.fileMeta.timePoints;
    n = P.samples(idx_sample).n_t_meshes;
    if ~isempty(n), timePoints = timePoints(1:n); end

    meshes = cell(length(timePoints), 1);
    radii = zeros(length(timePoints), nU);
    for tidx = 1:length(timePoints)
        tubi.setTime(timePoints(tidx));
        mesh = tubi.getCurrentSPCutMesh();
        mesh.v = tubi.xyz2APDV(mesh.v);                % to the anatomical frame
        v = reshape(mesh.v, nU, nV, 3);
        radii(tidx,:) = mean(sqrt(sum((v - mean(v, 2)).^2, 3)), 2);
        meshes{tidx} = struct('f', mesh.f, 'v', mesh.v);
    end
    save(fullfile(P.meshes, ['meshes_sample' num2str(idx_sample)]), 'meshes', 'radii');
    fprintf('sample %d: %d meshes\n', idx_sample, numel(meshes));
end
