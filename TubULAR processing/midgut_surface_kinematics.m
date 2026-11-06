%% MIDGUT SURFACE KINEMATICS (TubULAR) -- one copy of this script per sample
% Run after midgut_tubular_pipeline.m. Time-smoothed (s,phi) meshes and
% pullbacks -> PIV -> 3D velocities -> metric kinematics, then pathlines
% and strain from a reference time t0 set by the cross-sample time
% alignment (Supplement Secs. IV.3-IV.4, Fig. S7). "Methods" and
% "Supplement" refer to "Mechanical
% canalization of chiral gut morphogenesis". See README.md.


%% USER SETTINGS ==========================================================
clear; close all;

% Code
tubularDir   = '/path/to/tubular' ;   % same TubULAR as for midgut_tubular_pipeline.m
gptoolboxDir = '' ;                   % separate gptoolbox install (optional; TubULAR bundles one in external/)

% Data: the same sample as in midgut_tubular_pipeline.m
dataDir    = '/path/to/sample/deconvolved' ;           % holds masterSettings.mat
dir16bit   = fullfile(dataDir, 'deconvolved_16bit') ;  % 16-bit data directory = TubULAR project directory (holds xp.mat)
shiftsFile = '/path/to/timeshifts.mat' ;               % time shifts from the cross-sample alignment, read by alignedT0

overwrite = false ; % true: recompute steps even if their output already exists

addpath(genpath(tubularDir))
if ~isempty(gptoolboxDir), addpath(genpath(gptoolboxDir)) ; end


%% TubULAR INSTANCE =======================================================
load(fullfile(dataDir, 'masterSettings.mat'), 'masterSettings')
cd(dir16bit)
load(fullfile(dir16bit, 'xp.mat'), 'xp', 'opts')
tubi = TubULAR(xp, opts) ;
tubi.xp.expMeta.channelsUsed = masterSettings.channelsUsed + 1 ;  % TubULAR channel indices are 1-based


%% (Optional) Mesh evolution renderings
% tidx2do needs the visualizeMeshEvolution of our TubULAR version (see README).
options = struct('overwrite', overwrite, 'plot_growth', false, 'plot_meshOnly', false, ...
    'meshcolor', 0.7 * [1, 1, 1], 'planeColor', true, ...
    'tidx2do', 16) ;   % index into the timepoints to render
tubi.visualizeMeshEvolution(options)

options = struct('overwrite', overwrite, 'faceon', false, 'meshcolor', 0.7 * [1, 1, 1]) ;
tubi.visualizeMeshEvolution(options)


%% PART 1: TIME-SMOOTHED MESHES AND PULLBACKS =============================

%% Smooth the (s,phi) meshes in time
smOpts = struct('overwrite', overwrite, 'width', 2) ;  % kernel width, in timepoints
tubi.smoothDynamicSPhiMeshes(smOpts) ;

%% (Optional) Plot the smoothed meshes and their coordinate system
tubi.plotSPCutMeshSmRS(smOpts) ;
tubi.coordinateSystemDemo(struct('coordSys', 'spsm'))

%% Pullbacks of the time-smoothed meshes
% MIPs over 15 onion layers outward of the surface, 1.5 px apart, in all
% channels. As in our runs, this also regenerates the (s,phi) pullbacks
% made by midgut_tubular_pipeline.m with these settings; set
% generate_sphi = false to keep those instead.
for tt = tubi.xp.fileMeta.timePoints
    disp(['NOW PROCESSING TIME POINT ', num2str(tt)]);
    tubi.setTime(tt) ;
    tubi.getCurrentData() ;
    pbOptions = struct() ;
    pbOptions.overwrite     = true ;     % needed to replace the existing (s,phi) pullbacks
    pbOptions.numLayers     = [15 0] ;   % [outward, inward]
    pbOptions.layerSpacing  = 1.5 ;      % pixels between onion layers
    pbOptions.generate_spsm = true ;     % pullbacks of the time-smoothed meshes (PIV input)
    pbOptions.generate_sphi = true ;
    tubi.generateCurrentPullbacks([], [], [], pbOptions) ;
end

%% Tile the smoothed pullbacks in phi (images for PIV)
tubi.doubleCoverPullbackImages(struct('coordsys', 'spsm'))


%% PART 2: PIV AND EULERIAN KINEMATICS (Supplement Sec. IV.4) =============

%% PIV on the tiled pullbacks (PIVlab by default)
pivOpts = struct() ;
pivOpts.overwrite = overwrite ;
pivOpts.intArea1  = 256 ;   % interrogation window of pass 1, in pixels (TubULAR default 128)
pivOpts.intArea2  = 128 ;   % pass 2 (default 64); passes 3-4 keep the defaults (32, 16)
tubi.measurePIV2d(pivOpts) ;

%% 3D velocities
tubi.measurePIV3d(struct('overwrite', overwrite, 'show_v3d_on_data', false)) ;

%% Lagrangian (pathline) time-averaging of the velocities
tubi.timeAverageVelocities(struct('overwrite', overwrite))
options = struct() ;
options.plot_vxyz      = false ;
options.invertImage    = true ;
options.averagingStyle = 'Lagrangian' ;
tubi.plotTimeAvgVelocities(options)

%% Divergence and curl of the tangential flow (Helmholtz-Hodge)
options = struct() ;
options.overwrite      = overwrite ;
options.averagingStyle = 'Lagrangian' ;
options.lambda         = 0 ;
options.lambda_mesh    = 0 ;
tubi.helmholtzHodge(options) ;

%% Metric kinematics: div(v_parallel), 2Hv_n and area change rate (Eqs. S2-S3)
tubi.measureMetricKinematics(struct('overwrite', overwrite))
tubi.plotMetricKinematics(struct())


%% PART 3: PATHLINES FROM THE ALIGNED REFERENCE TIME ======================

%% Restrict the timepoints to those covered by PIV (N PIV frames <-> N+1 timepoints)
P = load(tubi.fileName.pivRaw.raw, 'u_filtered') ;
nKeep = min(numel(P.u_filtered) + 1, numel(tubi.xp.fileMeta.timePoints)) ;
tps = tubi.xp.fileMeta.timePoints(1:nKeep) ;
tubi.xp.fileMeta.timePoints = tps ;
if isstruct(tubi.xp), tubi.xp.tIdx = containers.Map(tps, 1:length(tps)) ; end
tubi.piv.smoothed = struct() ;
tubi.loadPIV() ;   % reload the PIV, sized to the trimmed list
fprintf('Using timepoints %d to %d (%d PIV frames)\n', tps(1), tps(end), numel(P.u_filtered)) ;

%% Reference time t0 from the cross-sample time alignment (Supplement Sec. IV.3)
% alignedT0 (from the time-alignment code) returns this sample's reference
% timepoint from the shifts in shiftsFile. Writing t0 to disk keeps
% TubULAR's internal t0set() calls from resetting it. Steps above use t0
% only for plots (labels and time axes).
t0 = alignedT0(tubi, shiftsFile) ;
write_txt_with_header(tubi.fileName.t0, t0, 't0, the reference timestamp') ;
tubi.t0set(t0) ;

%% Pathlines seeded at t0, and accumulated strain along them
tubi.generateRicciMeshTimePoint(t0, struct('radiusTolerance', 0.03)) ;  % conformal reference mesh at t0
tubi.measurePullbackPathlines(struct('overwrite', overwrite, 't0Pathlines', t0)) ;
tubi.measurePathlineStrain(struct('overwrite', overwrite, 't0Pathline', t0)) ;
tubi.plotPathlineStrain(struct('t0Pathline', t0)) ;

%% Velocities and metric kinematics along the pathlines
tubi.measurePathlineVelocities(struct('overwrite', overwrite))
tubi.plotPathlineVelocities(struct('gridTopology', 'triangulated'))
tubi.measurePathlineMetricKinematics(struct('overwrite', overwrite, 't0Pathline', t0))
tubi.plotPathlineMetricKinematics(struct('t0Pathline', t0))

%% Strain rates, and along the pathlines
tubi.measureStrainRate(struct('overwrite', overwrite))
tubi.plotStrainRate3DFiltered(struct())
tubi.plotStrainRate(struct('clim_trace', 0.05, 'clim_deviatoric', 0.05))
tubi.measurePathlineStrainRate(struct('overwrite', overwrite, 't0Pathline', t0, ...
    'overwriteImages', false, 'plot_dzdp', false))
tubi.plotPathlineStrainRate(struct('climit', 0.05, 'climitWide', 1.0))
