% extractNucleiStripes_noTracks.m
% Hand-tracing of the four longitudinal myocyte-nuclei stripes on the
% pullback images, and their embedding in 3D (Sec. IV.2 of the Supplement).
%
% 1. Draw: for every timepoint and each of the 4 stripes, draw the stripe on
%    the pullback image with the freehand tool; double-click to validate a
%    piece, press 'f' to go to the next frame, 'r' to redo the current frame.
%    Curves already on disk are skipped, so the drawing can be resumed.
%    Stripe order: 0 = red (VL), 1 = green (DL), 2 = yellow (DR), 3 = blue (VR).
% 2. Optional: an example pullback with the curves and a pullback movie.
% 3. Embed: resample each curve and push it to 3D on the (anatomical) cut
%    mesh, saving tubular_output/cellTracking/muscleStripes/muscleNucleiStripes.mat
%    (read by process_all_V6.m via load_myocytes_stripes.m).
%
% Requires TubULAR. Set idx_sample below.

clc; close all
clearvars -except tubi idx_sample
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'), here); P = config_paths();

if ~exist('idx_sample', 'var'), idx_sample = 1; end
if ~exist('tubi', 'var'), tubi = load_tubi(P, idx_sample); end

outdir = fullfile(tubi.dir.tracking, 'muscleStripes');
if isfield(P, 'stripes_dir') && ~isempty(P.stripes_dir)     % optional override (config_paths_local.m)
    outdir = fullfile(P.stripes_dir, sprintf('sample%d', idx_sample));
end
if ~exist(outdir, 'dir'), mkdir(outdir); end

isSmooth = false;          % draw on the smoothed-mesh pullbacks instead
example_pullback = true;   % figure of one pullback with the four curves
pullback_movie = true;     % movie of all pullbacks with the curves
preview = true;            % save 3D previews of the embedded curves
tp_example = 40;

if isSmooth
    imFileBase = tubi.fullFileBase.im_sp_sm;
else
    imFileBase = tubi.fullFileBase.im_sp;
end
% sample 8 was drawn on the pullbacks of its separate grid folder:
% imFileBase = fullfile(P.samples(8).tubular_output, 'gridCoords_nU0100_nV0100_forNico/PullbackImages_000step_sphi/Time_%06d_pbsp.tif');

ncurves = 4;
timePoints = tubi.xp.fileMeta.timePoints;
n = P.samples(idx_sample).n_t_meshes;
if ~isempty(n), timePoints = timePoints(1:n); end

%% 1. draw the stripes
global restartLoopFlag finishedFlag
restartLoopFlag = false;
finishedFlag = false;
fig = figure('KeyPressFcn', @myKeyPressCallback);
for curv = 0:ncurves-1
    curvfn = fullfile(outdir, sprintf('curvePositions_stripe%d.mat', curv));
    if exist(curvfn, 'file')
        load(curvfn, 'curvePositions')
        disp(['Curve ' num2str(curv) ' already exists. loaded from disk'])
    else
        curvePositions = cell(length(timePoints), 1);
    end

    tidx = 0;
    while tidx < length(timePoints)
        if ~restartLoopFlag, tidx = tidx + 1; end
        tp = timePoints(tidx);
        disp(['Timepoint ', num2str(tp)])
        if isempty(curvePositions{tidx}) || restartLoopFlag
            restartLoopFlag = false;
            clf
            im = imread(sprintf(imFileBase, tp));
            imshow(im(:,:,1));
            hold on;
            title('double click on ROI to validate')
            pos_rois = [];
            while ~finishedFlag
                h = drawfreehand;
                wait(h);
                pos_rois = [pos_rois; h.Position]; %#ok<AGROW>
            end
            curvePositions{tidx} = pos_rois;
            finishedFlag = false;
            disp('saving updates...')
            save(curvfn, 'curvePositions')
        end
    end
end

%% 2. example pullback and pullback movie with the curves
blend = @(im) cat(3, 0.95*255/255*im(:,:,1) + 0.9*102/255*im(:,:,2), ...
                     0.95*231/255*im(:,:,1) + 0.9*204/255*im(:,:,2), ...
                     0.95*146/255*im(:,:,1) + 0.9*255/255*im(:,:,2));
if example_pullback || pullback_movie
    C = cell(1, ncurves);
    for curv = 0:ncurves-1
        C{curv+1} = load(fullfile(outdir, sprintf('curvePositions_stripe%d.mat', curv)), 'curvePositions').curvePositions;
    end
end
if example_pullback
    figure
    imshow(blend(imread(sprintf(imFileBase, tp_example)))); hold on
    for curv = 1:ncurves
        XY = C{curv}{tp_example};
        plot(XY(:,1), XY(:,2), 'LineWidth', 8)
    end
end
if pullback_movie
    vidfile = VideoWriter(fullfile(outdir, sprintf('pullback_sample%d', idx_sample)));
    vidfile.FrameRate = 15;
    open(vidfile);
    figure
    for tp = 1:(length(timePoints)-1)
        clf
        imshow(blend(imread(sprintf(imFileBase, tp)))); hold on
        for curv = 1:ncurves
            XY = C{curv}{tp};
            plot(XY(:,1), XY(:,2), 'LineWidth', 8)
        end
        title(['frame tp=' num2str(tp)])
        img = imresize(frame2im(getframe(gcf)), [758, 1250]);
        writeVideo(vidfile, img);
    end
    close(vidfile);
end

%% 3. push the curves to 3D
close all
im = imread(sprintf(imFileBase, tubi.t0));
[~, ~, ~, xyzlims] = tubi.getXYZLims();
Xmax = size(im, 2);
nUi = 150;          % resampled points per curve
max_px = 1000;      % period used to unwrap the curves in Y
if preview
    for v = {'perspective2', 'dorsal'}
        if ~exist(fullfile(outdir, 'images', v{1}), 'dir'), mkdir(fullfile(outdir, 'images', v{1})); end
    end
end
dat = cell(length(timePoints), 1);
for tidx = 1:length(timePoints)
    tp = timePoints(tidx);
    tubi.setTime(tp);
    if isSmooth
        mesh = tubi.getCurrentSPCutMeshSmRS();
    else
        mesh = tubi.getCurrentSPCutMesh();
        mesh.v = tubi.xyz2APDV(mesh.v);                % to the anatomical frame
    end

    xyz = cell(4, 1);
    uvs = cell(4, 1);
    for curv = 0:ncurves-1
        load(fullfile(outdir, sprintf('curvePositions_stripe%d.mat', curv)), 'curvePositions')
        XY = curvePositions{tidx};
        [XX, uids] = unique(XY(:, 1));
        YY = XY(uids, 2);
        YY = max_px * unwrap(YY/max_px*2*pi) / (2*pi);

        % resample along the pullback's long axis and map to (u, v)
        Xi = linspace(1, Xmax, nUi);
        Yi = interp1(XX, YY, Xi, 'linear', 'extrap');
        uv = tubi.XY2uv(im, [Xi; Yi]', false, 1, 1);
        uv(:,2) = mod(uv(:,2), 1);
        push = uv(:, 1) == 0;                          % keep away from the boundary
        uv(push, 1) = eps;
        uv = uv(uv(:,1) < 1, :);

        uvs{curv+1} = uv;
        if isSmooth
            xyz{curv+1} = tubi.uv2APDV(uv, 'spsmrs', 1, 1);
        else
            xyz{curv+1} = tubi.uv2APDV(uv, 'sprs', 1, 1);
        end
        assert(~any(any(isnan(xyz{curv+1}))))
    end

    if preview
        close all;
        figure('position', [0 0 875 656])
        trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none', ...
            'facecolor', 0.8*[1,1,1], 'facealpha', 0.7);
        hold on;
        for k = 1:4
            plot3(xyz{k}(:,1), xyz{k}(:,2), xyz{k}(:,3), '-', 'linewidth', 4);
        end
        axis equal; axis off
        lighting gouraud; camlight
        xlim(xyzlims(1,:)); ylim(xyzlims(2,:)); zlim(xyzlims(3,:))
        set(gcf, 'color', 'w')
        view([48, 45]); drawnow
        title(['t=' num2str(tp*tubi.timeInterval) ' ' tubi.timeUnits])
        saveas(gcf, fullfile(outdir, 'images', 'perspective2', sprintf('muscleStripes_%06d.png', tp)))
        view(2)
        saveas(gcf, fullfile(outdir, 'images', 'dorsal', sprintf('muscleStripes_dorsal_%06d.png', tp)))
    end

    dat{tidx} = struct('mesh', mesh, ...
        'curv1', xyz{1}, 'curv2', xyz{2}, 'curv3', xyz{3}, 'curv4', xyz{4}, ...
        'uv1', uvs{1}, 'uv2', uvs{2}, 'uv3', uvs{3}, 'uv4', uvs{4});
    disp(['done with tp=' num2str(tp)])
end
save(fullfile(outdir, 'muscleNucleiStripes.mat'), 'dat')
