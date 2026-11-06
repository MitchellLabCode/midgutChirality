%% MIDGUT SURFACE PIPELINE (TubULAR) -- one copy of this script per sample
% Fused, deconvolved light-sheet stacks of the embryonic Drosophila midgut
% -> 16-bit stabilized stacks -> surface meshes (ilastik + level sets)
% -> APDV-aligned, parameterized cut meshes and pullbacks. Continue with
% midgut_tubular_kinematics.m for time-smoothed pullbacks, PIV and surface
% kinematics. Adapted from the TubULAR example pipelines (Mitchell & Cislo,
% Nat. Methods 2023). "Methods" and "Supplement" refer to "Mechanical canalization of chiral gut morphogenesis".
% See README.md for requirements.
%
% Run section by section. After restarting MATLAB, rerun USER
% SETTINGS and MASTER SETTINGS before continuing. Two steps happen in
% ilastik, where indicated.


%% USER SETTINGS ==========================================================
clear; close all;

% Code
tubularDir   = '/path/to/tubular' ;   % TubULAR, with its DECLab and TexturePatch submodules
gptoolboxDir = '' ;                   % separate gptoolbox install (optional; TubULAR bundles one in external/)

% Data
dataDir          = '/path/to/sample/deconvolved' ;           % fused, deconvolved 32-bit TIFFs (+ masterSettings.mat)
dir32bit         = dataDir ;
dir16bit         = fullfile(dataDir, 'deconvolved_16bit') ;  % 16-bit data directory = TubULAR project directory
dir16bit_prestab = fullfile(dir16bit, 'data_pre_stabilization') ;

% Sample-specific reference timepoints (timestamps, not indices)
trefAPDV    = 39 ;  % timepoint at which the dorsal landmark is trained in ilastik
trefEndcaps = 0 ;   % timepoint at which the dorsal-most endcap vertices are defined
t0_for_phi0 = 0 ;   % timepoint at which phi0 (the corners of the (s,phi) pullback maps) is defined
% The Lagrangian reference time t0 is set in midgut_tubular_kinematics.m.

overwrite = false ; % true: recompute steps even if their output already exists

addpath(genpath(tubularDir))
if ~isempty(gptoolboxDir), addpath(genpath(gptoolboxDir)) ; end
cd(dataDir)


%% MASTER SETTINGS ========================================================
% Saved to masterSettings.mat on first run and loaded from disk afterwards.
% After editing values here, delete that file (or run with overwrite = true).
masterSettingsFn = fullfile(dataDir, 'masterSettings.mat') ;
if ~exist(masterSettingsFn, 'file') || overwrite
    masterSettings = struct( ...
        'stackResolution', [0.195 0.195 0.195], ... % spaceUnits per pixel
        'nChannels', 3, ...
        'channelsUsed', [0 1 2], ...        % 0-based, as in the filenames
        'timePoints', 0:85, ...             % timepoints to pre-process
        'ssfactor', 5, ...                  % subsampling for ilastik and mesh extraction
        'flipy', true, ...                  % data volume is mirrored wrt the lab frame (sets handedness!)
        'timeInterval', 2, ...              % time between timepoints, in timeUnits
        'timeUnits', 'min', ...
        'spaceUnits', [char(956) 'm'], ...
        'scale', [3 0.1 0.7], ...           % per-channel intensity scale for 32 -> 16 bit
        'file32Base', 'TP%d_Ch%d_Ill0_Ang58,118,178,238,298,358.tif', ... % fused, deconvolved filenames
        'fn_prestab', 'Time_%06d_c%d.tif', ...  % 16-bit, before stabilization
        'fn', 'Time_%06d_c%d_stab', ...         % 16-bit, stabilized
        'swapZT', 0, ...
        'set_preilastikaxisorder', 'cxyz', ...  % axis order of the h5 files written for ilastik
        'nU', 100, ...                      % longitudinal samples of the surface grid
        'nV', 100) ;                        % circumferential samples of the surface grid
    save(masterSettingsFn, 'masterSettings')
else
    load(masterSettingsFn, 'masterSettings')
end

stackResolution = masterSettings.stackResolution ;
nChannels       = masterSettings.nChannels ;
channelsUsed    = masterSettings.channelsUsed ;
timePoints      = masterSettings.timePoints ;
ssfactor        = masterSettings.ssfactor ;
flipy           = masterSettings.flipy ;
timeInterval    = masterSettings.timeInterval ;
timeUnits       = masterSettings.timeUnits ;
spaceUnits      = masterSettings.spaceUnits ;
scale           = masterSettings.scale ;
file32Base      = masterSettings.file32Base ;
fn_prestab      = masterSettings.fn_prestab ;
fn              = masterSettings.fn ;
swapZT          = masterSettings.swapZT ;
set_preilastikaxisorder = masterSettings.set_preilastikaxisorder ;
nU              = masterSettings.nU ;
nV              = masterSettings.nV ;


%% PART 0: PRE-PROCESSING (Methods 5) =====================================
% 32-bit stacks -> 16-bit, jitter-corrected, multichannel TIFFs.

%% Convert 32-bit to 16-bit (existing outputs are skipped)
convert32to16bit(timePoints, scale, dir32bit, dir16bit_prestab, ...
    file32Base, fn_prestab, channelsUsed)

%% Maximum intensity projections (used to measure jitter)
mipDir = fullfile(dir16bit_prestab, 'mips') ;
mipOpts = struct() ;
mipOpts.overwrite_mips = overwrite ;
mipOpts.scale = -1 ;                % do not rescale intensities
mipOpts.channels = channelsUsed ;
makeMips(timePoints, dir16bit_prestab, fn_prestab, mipDir, mipOpts)

%% Collate channels into one multichannel TIFF per timepoint (pre-stabilization)
collateColors(fullfile(dir16bit_prestab, fn_prestab), ...
    fullfile(dir16bit_prestab, 'Time_%06d.tif'), timePoints, channelsUsed) ;

%% Remove translational jitter
% Stage/sample motion is measured in XYZ by phase correlation of the MIPs
% and subtracted; rotations are not corrected. The shifts are computed on
% the first channel in channelsUsed and reused for the other channels.
fileNameIn      = fullfile(dir16bit_prestab, fn_prestab) ;
fileNameOut     = fullfile(dir16bit, [fn '.tif']) ;
mipoutdir       = fullfile(mipDir, 'mips_stab') ;   % stabilized MIPs
mips_stab_check = fullfile(mipDir, 'stab_check') ;  % overlays for checking the stabilization
stabOpts = struct() ;
stabOpts.overwrite_tiffs = overwrite ;
stabOpts.im_intensity    = 1 ;      % intensity scale of the output MIPs
stabOpts.imref_intensity = 1 ;      % intensity scale of the reference MIP
for ch = channelsUsed
    stabOpts.stabChannel    = ch ;
    stabOpts.overwrite_mips = overwrite && ch == channelsUsed(1) ;
    stabilizeImages(fileNameIn, fileNameOut, [fn '.png'], 'uint16', ...
        timePoints, timePoints, timePoints(1), ...
        mipDir, mipoutdir, mips_stab_check, stabOpts)
end

%% Collate stabilized channels into multichannel TIFFs (TubULAR input)
collateColors(fullfile(dir16bit, [fn '.tif']), fullfile(dir16bit, 'Time_%06d.tif'), ...
    timePoints, channelsUsed) ;
% The single-channel *_stab.tif files are not needed after this step.


%% PART 1: TubULAR PROJECT AND SURFACE SEGMENTATION (Supplement Sec. IV.1)

%% Define the experiment (saved to xp.mat in the 16-bit data directory)
cd(dir16bit)
projectDir = dir16bit ;
fn = 'Time_%06d' ;                                            % collated, stabilized stacks
mslsDir = [fullfile(projectDir, 'tubular_output') filesep] ;  % TubULAR output directory

xpFn = fullfile(projectDir, 'xp.mat') ;
if ~exist(xpFn, 'file') || overwrite
    fileMeta                 = struct() ;
    fileMeta.dataDir         = dir16bit ;
    fileMeta.filenameFormat  = [fn '.tif'] ;
    fileMeta.nChannels       = nChannels ;
    fileMeta.timePoints      = timePoints ;   % may be restricted, e.g. to drop late timepoints
    fileMeta.stackResolution = stackResolution ;
    fileMeta.swapZT          = swapZT ;

    expMeta                  = struct() ;
    expMeta.channelsUsed     = channelsUsed ;
    expMeta.channelColor     = 1 ;
    expMeta.description      = 'Drosophila midgut' ;
    expMeta.dynamicSurface   = 1 ;
    expMeta.jitterCorrection = 0 ;            % jitter was removed in Part 0
    expMeta.fitTime          = fileMeta.timePoints(1) ;

    % Surface detection: level sets (MATLAB activecontour, Chan-Vese) on the
    % ilastik probabilities. Parameters marked * are tuned per sample.
    detectOptions = struct( ...
        'channel', 1, ...
        'foreGroundChannel', 2, ...       % ilastik probability channel of the gut label
        'ssfactor', ssfactor, ...
        'niter', 1000, ...                % * iterations per timepoint
        'niter0', 1000, ...               % * iterations for the first timepoint
        'pre_pressure', 20, ...           % * applied to the initial guess before evolving
        'pre_tension', 0, ...
        'pressure', 0, ...                % *
        'tension', 0.5, ...               % *
        'post_pressure', 2, ...
        'post_tension', 3, ...
        'exit_thres', 1e-7, ...
        'radius_guess', 15, ...           % * initial guess radius, in subsampled pixels
        'center_guess', 'click', ...      % * click the initial guess center, or 'x,y,z' in subsampled pixels
        'init_ls_fn', 'msls_initguess.mat', ...
        'timepoint', timePoints(1), ...
        'zdim', 2, ...
        'fileName', fn, ...
        'mslsDir', mslsDir, ...
        'ofn_ls', 'msls_', ...
        'ofn_ply', 'mesh_ms_', ...
        'ofn_smoothply', 'mesh_', ...
        'run_full_dataset', projectDir, ...
        'dset_name', 'exported_data', ...
        'save', true, ...
        'plot_mesh3d', false, ...
        'previewIsosurface', false, ...
        'dtype', 'h5', ...
        'mask', 'none', ...
        'mesh_from_pointcloud', false, ...
        'prob_searchstr', '_stab_Probabilities.h5', ...
        'preilastikaxisorder', set_preilastikaxisorder, ...
        'ilastikaxisorder', 'xyzc', ...
        'physicalaxisorder', 'xyzc', ...
        'include_boundary_faces', true, ...
        'smooth_with_matlab', 0.01, ...   % >0: MATLAB Laplacian smoothing with this lambda; <0: MeshLab (then also set 'mlxprogram')
        'useAxesPixels4Permutation', true) ;  % permute the ilastik output to xyzc by matching axis lengths (see below)
    save(fullfile(projectDir, 'msls_detectOpts.mat'), 'detectOptions')

    xp = struct('fileMeta', fileMeta, 'expMeta', expMeta, 'detectOptions', detectOptions) ;

    opts = struct() ;
    opts.meshDir      = mslsDir ;
    opts.flipy        = flipy ;
    opts.timeInterval = timeInterval ;
    opts.timeUnits    = timeUnits ;
    opts.spaceUnits   = spaceUnits ;
    opts.nU           = nU ;
    opts.nV           = nV ;
    opts.normalShift  = 0 ;           % dilation along surface normals for texture mapping
    opts.a_fixed      = 2.0 ;         % aspect ratio of pullback images
    opts.adjustlow    = 1.00 ;        % lower intensity percentile (clip)
    opts.adjusthigh   = 99.9 ;        % upper intensity percentile (clip)
    opts.phiMethod    = 'texture' ;   % surface-Lagrangian tracking of phi (TubULAR default: 'curves3d')
    opts.lambda_mesh  = 0.00 ;        % mesh smoothing before DEC measurements
    opts.lambda       = 0.0 ;         % smoothing of measured fields on the surface
    opts.lambda_err   = 0.0 ;
    opts.zwidth       = 1 ;
    opts.nmodes       = 7 ;
    save(xpFn, 'xp', 'opts')
else
    load(xpFn, 'xp', 'opts')
end

%% TubULAR instance
tubi = TubULAR(xp, opts) ;
tubi.xp.expMeta.channelsUsed = channelsUsed + 1 ;  % TubULAR channel indices are 1-based
tubi.plotting.preview = false ;
tps = tubi.xp.fileMeta.timePoints(:)' ;

%% Write subsampled h5 stacks for ilastik
tubi.prepareIlastik() ;

%% Record the axis lengths of the ilastik input
% The h5 files are written in cxyz order so that ilastik's three orthogonal
% views display properly. getMeshes then permutes the exported
% probabilities back to xyzc by matching axis lengths, so x, y and z must
% have different lengths after subsampling.
h5fn = fullfile(tubi.dir.data, [sprintf(fn, tps(1)) '.h5']) ;
[cnum, xpix, ypix, zpix] = size(h5read(h5fn, '/inputData')) ;
if numel(unique([xpix, ypix, zpix])) < 3
    error('x, y and z have equal lengths (%d, %d, %d) after subsampling; pad the data so that they differ.', ...
        xpix, ypix, zpix)
end
tubi.axesPixels = struct('xpix', xpix, 'ypix', ypix, 'zpix', zpix, 'cnum', cnum) ;
tubi.axes_pixels_list = [xpix, ypix, zpix, cnum] ;

%% ilastik: segment the gut (outside MATLAB)
% Train an Autocontext (2-stage) project on the h5 files with four labels:
% the gut (including the muscle layers and enclosed yolk), background, and
% anterior and posterior endcaps that close the gut at its ends for the
% level-set segmentation. Export the probabilities of every timepoint as
% Time_XXXXXX_Probabilities.h5 next to the h5 inputs (training on a subset
% of timepoints is usually enough). detectOptions.foreGroundChannel must
% point to the gut label.

%% Extract meshes by level sets on the ilastik probabilities
% To re-tune the detection parameters after xp.mat exists, edit
% tubi.xp.detectOptions here before calling getMeshes.
tubi.xp.detectOptions.timepoints = fliplr(tps) ;   % last timepoint first
tubi.getMeshes(overwrite)


%% PART 2: APDV FRAME AND SURFACE PARAMETERIZATION (Methods 6; Supplement Sec. IV.2)

%% ilastik: APDV landmarks (outside MATLAB)
% In a separate ilastik project, train blobs on the anterior end (junction
% with the foregut), the posterior end, a dorsal point and background, each
% in its own label channel; A and P at every timepoint, D only at trefAPDV.
% Export as Time_XXXXXX_Probabilities_apcenterline.h5. The label order must
% match anteriorChannel and posteriorChannel below.

%% Global APDV frame
alignAPDVOpts = struct() ;
alignAPDVOpts.overwrite = overwrite ;
alignAPDVOpts.use_iLastik = true ;
alignAPDVOpts.ilastikOutputAxisOrder = 'cxyz' ;
alignAPDVOpts.tref = trefAPDV ;
tubi.computeAPDVCoords(alignAPDVOpts) ;

%% Anterior and posterior points at every timepoint
apdOpts = struct() ;
apdOpts.overwrite = overwrite ;
apdOpts.use_iLastik = true ;
apdOpts.ilastikOutputAxisOrder = 'cxyz' ;
apdOpts.anteriorChannel  = 2 ;
apdOpts.posteriorChannel = 3 ;
apdOpts.thres = 0.99 ;
apdOpts.preview = true ;
apdOpts.preview_com = false ;
apdOpts.preview_probabilities = false ;
[apts_sm, ppts_sm] = tubi.computeAPDpoints(apdOpts) ;

%% Align the meshes in the APDV frame
% tubi.flipy mirrors the data into the lab frame from here on, so it sets
% the handedness of every downstream result: verify it on a wild-type sample.
alignAPDVOpts.overwrite_ims = overwrite ;
tubi.alignMeshesAPDV(alignAPDVOpts) ;

%% (Optional, slow) Texture-mapped 3D renderings of the aligned meshes
metadat = struct() ;
metadat.timePoints         = tps([1 end]) ;  % subset of timepoints to render
metadat.reorient_faces     = false ;
metadat.overwrite          = overwrite ;
metadat.normal_shift       = tubi.normalShift ;  % push along normals, in data-space pixels
metadat.texture_axis_order = [1 2 3] ;           % axis permutation between surface and data, if any
metadat.smoothing_lambda   = 0.001 ;
metadat.perspective_angle  = [-195, 19] ;
metadat.plot_dorsal   = false ;
metadat.plot_ventral  = false ;
metadat.plot_lateral1 = false ;
metadat.plot_lateral2 = false ;
texOpts = struct() ;
texOpts.PSize        = 5 ;          % texture refinement on each face (higher is slower)
texOpts.numLayers    = [12, 12] ;   % layers to MIP over, [outward, inward]
texOpts.layerSpacing = 1.2 ;        % spacing between layers, in pixels
texOpts.smoothIter   = 0 ;
texOpts.isRGB        = false ;
texOpts.isFalseColor = true ;
texOpts.falseColors  = [1 0.78 0.39; 0 0 0; 0.39 0.706 1] ;  % one row per channel
tubi.plotSeriesOnSurfaceTexturePatch(metadat, texOpts)

%% Fast-marching centerlines
% Only used to constrain the cut path for the mapping to the plane, not for
% geometric measurements.
cntrlineOpts = struct() ;
cntrlineOpts.overwrite      = overwrite ;
cntrlineOpts.overwrite_ims  = overwrite ;
cntrlineOpts.weight         = 1e-11 ;  % speed-up weight; larger is faster but less precise
cntrlineOpts.exponent       = 1 ;      % scaling of the distance transform for speed through voxels
cntrlineOpts.res            = 1 ;      % resolution of the distance-transform grid
cntrlineOpts.preview        = true ;
cntrlineOpts.reorient_faces = false ;
cntrlineOpts.dilation       = 0 ;      % voxels to dilate the segmentation before computing paths
tubi.generateFastMarchingCenterlines(cntrlineOpts)

%% Flag anomalous centerlines
tubi.cleanFastMarchingCenterlines(struct('ssr_thres', 2000)) ;  % sum of squared residuals threshold (spaceUnits)

%% Cylinder meshes: remove the endcaps (topological sphere -> cylinder)
endcapOpts = struct('adist_thres', 10, ...   % distance threshold for cutting off the anterior, in pixels
                    'pdist_thres', 10, ...   % distance threshold for cutting off the posterior, in pixels
                    'tref', trefEndcaps) ;
if exist(tubi.fileName.endcapOptions, 'file') && ~overwrite
    tubi.loadEndcapOptions() ;
    endcapOpts = tubi.endcapOptions ;
else
    tubi.setEndcapOptions(endcapOpts) ;
    tubi.saveEndcapOptions() ;
end
methodOpts = struct('overwrite', overwrite, 'save_figs', true, 'preview', false) ;
tubi.sliceMeshEndcaps(endcapOpts, methodOpts) ;

%% Clean the cylinder meshes (remove 'ears' at the endcaps)
tubi.cleanCylMeshes(struct('overwrite', overwrite))

%% Cut meshes, (s,phi) meshes and pullbacks
% phi0 of each timepoint is matched to the (s,phi) pullback of its
% neighbor, so timepoints are processed outward from t0_for_phi0.
for tt = [tps(tps >= t0_for_phi0), fliplr(tps(tps < t0_for_phi0))]
    disp(['NOW PROCESSING TIME POINT ', num2str(tt)]);
    tubi.setTime(tt) ;

    % Cut the cylinder along a dorsal path joining its two ends
    cutMeshfn = sprintf(tubi.fullFileBase.cutMesh, tt) ;
    cutPathfn = sprintf(tubi.fullFileBase.cutPath, tt) ;
    if ~exist(cutMeshfn, 'file') || ~exist(cutPathfn, 'file') || overwrite
        tubi.generateCurrentCutMesh(struct())
        tubi.plotCutPath(tubi.currentMesh.cutMesh, tubi.currentMesh.cutPath)
        compute_pullback = true ;
    else
        compute_pullback = false ;
    end
    tubi.getCurrentUVCutMesh() ;

    % (s,phi) parameterization
    spcutMeshOptions = struct() ;
    spcutMeshOptions.t0_for_phi0     = t0_for_phi0 ;
    spcutMeshOptions.save_phi0patch  = false ;
    spcutMeshOptions.iterative_phi0  = false ;
    spcutMeshOptions.smoothingMethod = 'savgol' ;
    spcutMeshOptions.overwrite       = overwrite ;
    tubi.generateCurrentSPCutMesh([], spcutMeshOptions) ;

    % Pullback images: MIPs over onion layers around the surface. The
    % (s,phi) pullbacks are regenerated later by midgut_tubular_kinematics.m.
    if compute_pullback || ~exist(sprintf(tubi.fullFileBase.im_sp, tt), 'file') ...
            || ~exist(sprintf(tubi.fullFileBase.im_uv, tt), 'file')
        pbOptions = struct() ;
        pbOptions.overwrite      = overwrite ;
        pbOptions.generate_uv    = true ;    % conformal
        pbOptions.generate_sphi  = true ;    % (s,phi)
        pbOptions.generate_spsm  = false ;   % made in midgut_tubular_kinematics.m
        pbOptions.generate_ricci = false ;
        pbOptions.numLayers      = [20, 20] ;  % [outward, inward]
        pbOptions.layerSpacing   = 5 ;         % pixels between onion layers
        pbOptions.smoothIter     = 0 ;
        pbOptions.channels       = [1, 3] ;    % TubULAR channels: CVM nuclei and endoderm nuclei
        pbOptions.save_as_stack  = false ;
        tubi.generateCurrentPullbacks([], [], [], pbOptions) ;
    end
end

%% (Optional) Inspect the coordinate systems
tubi.setTime(tps(end)) ;
tubi.coordinateSystemDemo(struct('coordSys', 'uv'))

%% Export the (s,phi) cut meshes as OBJ files (e.g. for rendering in Blender)
% Numbered by timepoint index: obj_0.obj, obj_1.obj, ...
for tidx = 1:numel(tps)
    tubi.setTime(tps(tidx)) ;
    mesh = tubi.loadCurrentSPCutMesh() ;
    meshout = struct() ;
    meshout.f  = mesh.f ;
    meshout.v  = mesh.v ;
    meshout.vn = mesh.vn ;
    meshout.u  = mesh.sphi ;
    meshout.u(:, 1) = meshout.u(:, 1) / max(meshout.u(:, 1)) ;  % normalize s to [0, 1]
    exportOBJ(meshout, fullfile(tubi.dir.spcutMesh, sprintf('obj_%d.obj', tidx - 1)))
end


% Continue with midgut_tubular_kinematics.m (time-smoothed pullbacks, PIV,
% surface kinematics and pathlines).
