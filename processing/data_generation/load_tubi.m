function tubi = load_tubi(P, idx_sample)
% load_tubi  TubULAR object of one embryo, with its paths updated to the
% location of the data on this machine.
%
%   P = config_paths();  tubi = load_tubi(P, idx_sample);
%
% Requires TubULAR (P.tubular). Applies the per-sample overrides of
% config_paths.m (sample 8's cut meshes live in a separate grid folder).

addpath(genpath(P.tubular))
s = P.samples(idx_sample);
fprintf('loading xp struct of sample %d from %s\n', idx_sample, s.raw_dir);
load(fullfile(s.raw_dir, 'xp.mat'), 'xp', 'opts')

% the paths stored in xp.mat point to where the data were processed
xp.detectOptions.mslsDir = s.tubular_output;
xp.fileMeta.dataDir = s.raw_dir;
opts.meshDir = s.tubular_output;
tubi = TubULAR(xp, opts);

if ~isempty(s.spcutmesh_override)
    tubi.fullFileBase.spcutMesh = s.spcutmesh_override;
end
end
