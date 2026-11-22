%% SIFig_AreaConservation.m
%
% The four panels of SI Fig 6, for WT and bynMyo1c:
%   - each genotype's surface area and volume against aligned time
%     (ensemble mean with a shaded band across embryos)
%   - each genotype's areal rate of change, 1/2 Tr[g^-1 gdot] =
%     div(v_par) - 2Hv_n, against dilatational flow div(v_par) (mean of
%     the embryos' bin means, with SEM error bars over a band of +/- 1 SD)
% Every figure comes in three versions with the same size and axes
% placement, so they overlay exactly: labeled (axis labels, tick numbers,
% title, legend), ticks (tick numbers only) and bare (no text at all).
%   <outdir>/ensemble_<genotype>_area_volume_<style>.pdf/.png
%   <outdir>/arealRate_vs_divv_<genotype>_average_<style>.pdf/.png
% One TubULARCollection per genotype builds the TubULAR instances. Each
% specimen gets a folder, <outdir>/<genotype>/<YYYYMMDD>_<genotype>/.
%
% Time alignment
%   Each specimen's lag onto the WT reference's timeline is its
%   relative_shifts entry in shiftsFile, plus cross_genotype_shift for
%   bynMyo1c, so WT and bynMyo1c share one timeline:
%     aligned timepoint = timepoint + shiftSign * lag
%
% Surface area and volume
%   Each specimen's surface area and volume are divided by their value at
%   its aligned t0. xp.mat can list timepoints beyond the last aligned
%   mesh, so each specimen's timepoint list is first trimmed to end there.
%   The mean, SD and SEM at each time are over the specimens that have
%   data there. 
%
% Areal rate of change vs dilatational flow
%   From the measurements folder measureMetricKinematics wrote
%   (divv_vertices_*.mat, gdot_vertices_*.mat): every vertex in the middle
%   of the AP axis at every timepoint inside rateWindow (by default the
%   aligned window that every specimen covers) is binned by its div(v_par),
%   and its areal rate of change is averaged per bin, as in TubULAR's (Mitchell and Cislo, 2022, Nat. Methods)
%   Extended Data Fig. 4d. These are the raw per-vertex values, with no
%   quadrant averaging. Hoops whose areal rate is far off are left out
%   first. Each genotype's point is the mean of its specimens' bin means,
%   with SEM and SD across specimens. The bins span the central 99% of
%   div(v_par) at every timepoint on disk (rateBinRange = 'all') or inside
%   rateWindow ('window'); the console prints both ranges. Rates are in
%   1/min, converted with minutesPerTimepoint.
%
% Prerequisites (See TubULAR Processing subfolder)
%   - alignMeshesAPDV has been run for every dataset (aligned meshes exist
%     under <dataDir>/tubular_output).
%   - segment_strains_updated has been run for every dataset: the
%     TubULARCollection constructor requires DVLRstrain.mat, even though
%     these figures do not use it.
%   - @TubULAR/measureSurfaceAreaVolume.m is the patched version that
%     accepts options.outdir (the unpatched one writes into
%     <dataDir>/tubular_output instead; the script warns you).
%   - shiftsFile holds relative_shifts, cross_genotype_shift, dataDirs,
%     genotypes and <genotype>_indices, and every dataset it lists has its
%     xp.mat, since t = 0 depends on all of them.
%   - measureMetricKinematics has been run for every dataset.

addpath(genpath('.../gptoolbox'))
addpath(genpath('.../tubular'))

%% -------------------------------------------------------------------------
%  Settings
% --------------------------------------------------------------------------
outdir = 'out_dir';

overwrite = false;   % true = measure surface area and volume again, even if cached

% Time alignment: the strain-rate script's lags, with alignedT0.m's t = 0
shiftsFile          = 'shiftfile.mat';
shiftSign           = +1;   % aligned timepoint = timepoint + shiftSign * lag (check_shift_sign.m)
minutesPerTimepoint = 2;    % true time between frames [min], for every time axis; saved
                            % rates are per tubi.timeInterval and get rescaled to 1/min with this

% Areal rate of change vs dilatational flow, from the measurements folder
% measureMetricKinematics wrote in <dataDir>/tubular_output/
% gridCoords_nU0100_nV0100/metricKinematics/<mkRunFolder>/measurements
% ('' = the run matching the smoothing settings in each dataset's xp.mat)
mkRunFolder        = 'lambda0p000_lmesh0p000_lerr0p000_modes07w01';
useFiltered        = true;       % true: divv_filt & gdot_filt; false: unfiltered divv & gdot
rateWindow         = 'common';   % aligned times used: 'common' = the window every specimen
                                 % covers, 'all' = every timepoint on disk, or [t1 t2] in min
rateBinRange       = 'all';      % bins span the central 99% of div(v_par): 'all' = at every
                                 % timepoint on disk, as the earlier script did; 'window' = only
                                 % inside rateWindow (ignored when rateBinEdges is given)
apTrim             = 0.1;        % ignore this fraction of the AP axis at each end
zThresh            = 30;         % leave out hoops whose areal rate is more than zThresh robust
minDeviation       = 0.1;        % SDs and minDeviation [1/min] from typical (Inf = keep all)
rateBinEdges       = [];         % div(v_par) bin edges [1/min], or [] for nRateBins equal
nRateBins          = 25;         % bins spanning the central 99% of all specimens' div(v_par)
minPointsPerBin    = 20;         % a specimen's bin needs this many vertices
minSpecimensPerBin = 3;          % a genotype's point needs this many specimens
rateXLim           = [];         % axes shared by both areal-rate figures [1/min]:
rateYLim           = [];         % [] = automatic (the values used are printed)

% Surface area and volume
ensembleEndTime     = 110;      % latest end of both panels [min] ([] = no limit). They end
                                % together: at the last time every genotype still has
                                % minSpecimensForMean specimens, or here if that is earlier
minSpecimensForMean = 3;        % fewest specimens behind a mean; SD/SEM are over the specimens
                                % that have data at each time ([] = every specimen)
ensembleBand        = 'SD';     % shaded band across specimens: 'SD' or 'SEM'
ensembleYLim        = [];       % y axis shared by both figures: [] = automatic

if isempty(outdir)
    error('Set outdir at the top of the script before running.');
end
if ~exist(shiftsFile, 'file')
    error('shiftsFile not found: %s', shiftsFile);
end
if ischar(rateWindow)
    rateWindow = lower(rateWindow);
end
if ~((ischar(rateWindow) && any(strcmp(rateWindow, {'common', 'all'}))) || ...
        (isnumeric(rateWindow) && numel(rateWindow) == 2))
    error('rateWindow must be ''common'', ''all'' or [t1 t2] in min.');
end
rateBinRange = lower(rateBinRange);
if ~any(strcmp(rateBinRange, {'window', 'all'}))
    error('rateBinRange must be ''window'' or ''all''.');
end
if ~isempty(ensembleEndTime) && ~(isscalar(ensembleEndTime) && ensembleEndTime >= 0)
    error('ensembleEndTime must be [] or a time >= 0 in min.');
end
ensembleBand = upper(ensembleBand);
if ~any(strcmp(ensembleBand, {'SD', 'SEM'}))
    error('ensembleBand must be ''SD'' or ''SEM''.');
end
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

%% -------------------------------------------------------------------------
%  Dataset paths
% --------------------------------------------------------------------------
wtDirs = { ...
    'location_of_Tubular_Output_for_WT_sample_1', ...
    'location_of_Tubular_Output_for_WT_sample_2', ...
    };

bynDirs = { ...
     'location_of_Tubular_Output_for_bynMyo1C_sample_1', ...
    'location_of_Tubular_Output_for_bynMyo1C_sample_2', ...
    };

% 'name' is both the genotype folder under outdir and the suffix of each
% specimen folder (e.g. 20240528_bynmyo1c). 'title' is the figure title
% (TeX markup: \it ... \rm sets italics).
genotypes = struct( ...
    'name',  {'WT',        'bynmyo1c'}, ...
    'title', {'Wild-type', '\itbyn\rm > Myo1C'}, ...
    'dirs',  {wtDirs,      bynDirs});

% Settings for reading the metric kinematics
R = struct();
R.apTrim              = apTrim;
R.zThresh             = zThresh;
R.minDeviation        = minDeviation;
R.minutesPerTimepoint = minutesPerTimepoint;
if useFiltered
    R.dVar = 'divv_filt';
    R.gVar = 'gdot_filt';
else
    R.dVar = 'divv';
    R.gVar = 'gdot';
end
styles = {'labeled', 'ticks', 'bare'};   % every figure is saved in these three versions

%% -------------------------------------------------------------------------
%  Time alignment: every dataset's lag, and t = 0 
% --------------------------------------------------------------------------
TL = alignedTimeline(load(shiftsFile), shiftSign);
fprintf(['\nt = 0 is aligned timepoint %g, the first at which all %d datasets in the ' ...
    'shifts file have data\n(the last to start is dataset %d, %s).\n'], ...
    TL.alignedTimepointAtZero, numel(TL.lags), TL.lastToStart, TL.genotypes{TL.lastToStart});

%% -------------------------------------------------------------------------
%  Every specimen: surface area and volume, and the metric kinematics
% --------------------------------------------------------------------------
failed    = cell(0, 2);                  % rows of {specimen name, error message}
allCurves = cell(numel(genotypes), 1);   % per genotype, one curve per specimen
allRates  = cell(numel(genotypes), 1);   % per genotype, one set of vertex rates per specimen

for gg = 1:numel(genotypes)
    gName   = genotypes(gg).name;
    gDirs   = genotypes(gg).dirs;
    nSpec   = numel(gDirs);
    gOutdir = fullfile(outdir, gName);

    fprintf('\n========== %s: %d specimens ==========\n', gName, nSpec);

    % Folder name per specimen: <date>_<genotype>. makeUniqueStrings adds a
    % numeric suffix if two specimens of this genotype share a date.
    specNames = cell(nSpec, 1);
    for ii = 1:nSpec
        specNames{ii} = [specimenDate(gDirs{ii}) '_' gName];
    end
    specNames = matlab.lang.makeUniqueStrings(specNames);
    allCurves{gg} = cell(1, nSpec);   % stays empty for specimens that fail
    allRates{gg}  = cell(1, nSpec);

    % The constructor builds one TubULAR instance per dataset in col.tubis
    % (its own outdir is only written by measureRelativeTimeShifts)
    col = TubULARCollection(gDirs, repmat({gName}, nSpec, 1), gOutdir);

    for ii = 1:nSpec
        specOutdir = fullfile(gOutdir, specNames{ii});
        if ~exist(specOutdir, 'dir')
            mkdir(specOutdir);
        end
        fprintf('\n[%s %d/%d] %s\n  data: %s\n  out:  %s\n', ...
            gName, ii, nSpec, specNames{ii}, gDirs{ii}, specOutdir);
        tubi = col.tubis{ii};

        % Where this specimen sits on the aligned timeline
        try
            kData = findDataset(TL, gDirs{ii});
        catch ME
            fprintf(2, '  FAILED: %s\n', ME.message);
            failed(end+1, :) = {specNames{ii}, ['time alignment: ' ME.message]}; %#ok<AGROW>
            continue   % without a t0, nothing of this specimen can go on the timeline
        end
        t0 = datasetT0(TL, kData, shiftSign);
        fprintf('  lag %+g timepoints (dataset %d in the shifts file, %s), aligned t0 = %d\n', ...
            TL.lags(kData), kData, TL.genotypes{kData}, t0);
        xpMinutes = tubi.timeInterval * minutesPerUnit(tubi.timeUnits);
        if abs(xpMinutes - minutesPerTimepoint) > 1e-6 * minutesPerTimepoint
            fprintf(['  note: xp.mat gives %g %s per timepoint; time axes and rates here ' ...
                'use minutesPerTimepoint = %g min\n'], tubi.timeInterval, tubi.timeUnits, ...
                minutesPerTimepoint);
        end

        % Surface area and volume (measureSurfaceAreaVolume caches them)
        opts = struct();
        opts.overwrite = overwrite;
        opts.preview   = false;      % skip per-timepoint trisurf rendering
        opts.outdir    = specOutdir;
        try
            trimToAlignedMeshes(tubi);   % end timepoints at the last aligned mesh
            tps = tubi.xp.fileMeta.timePoints;
            if ~ismember(t0, tps)
                error(['The aligned t0 = %d is outside the timepoints with aligned meshes ' ...
                    '(%d to %d).'], t0, tps(1), tps(end));
            end
            tubi.t0 = t0;   % in memory only: measureSurfaceAreaVolume normalizes at tubi.t0
            tubi.measureSurfaceAreaVolume(opts);
            if ~exist(fullfile(specOutdir, 'surfacearea_volume_stab.mat'), 'file')
                error(['No surfacearea_volume_stab.mat in %s. MATLAB is probably using an ' ...
                    'unpatched measureSurfaceAreaVolume; check with ' ...
                    '"which measureSurfaceAreaVolume -all".'], specOutdir);
            end
            allCurves{gg}{ii} = loadCurve(tubi, specOutdir, specNames{ii}, gDirs{ii}, ...
                TL.lags(kData), t0, minutesPerTimepoint);
        catch ME
            fprintf(2, '  FAILED: %s\n', ME.message);
            failed(end+1, :) = {specNames{ii}, ME.message}; %#ok<AGROW>
        end
        close all   % measureSurfaceAreaVolume leaves its figure open

        % Vertex values for the areal rate of change vs dilatational flow,
        % each with its aligned time
        try
            allRates{gg}{ii} = loadRates(tubi, mkRunFolder, specNames{ii}, specOutdir, t0, R);
        catch ME
            fprintf(2, '  areal rate vs flow FAILED: %s\n', ME.message);
            failed(end+1, :) = {specNames{ii}, ['areal rate vs flow: ' ME.message]}; %#ok<AGROW>
        end
    end
end

%% -------------------------------------------------------------------------
%  Areal rate of change vs dilatational flow: each genotype's average over
%  the aligned times inside rateWindow
% --------------------------------------------------------------------------
rateWin = [];            % the aligned time window, set below
rates = [allRates{:}];   % every specimen slot, all genotypes
rates = [rates{:}];      % struct array of the specimens that loaded
xAll  = {};              % each specimen's div(v_par) at every timepoint on disk
xWin  = {};              % the same inside the window
if ~isempty(rates)
    % Aligned time window, as for the strain-rate script's figures 1-2
    if ischar(rateWindow) && strcmp(rateWindow, 'common')
        ranges  = vertcat(rates.tauRange);
        rateWin = [max(ranges(:, 1)), min(ranges(:, 2))];
        if rateWin(1) > rateWin(2)
            warning(['The specimens'' metric kinematics never all overlap in aligned time; ' ...
                'the areal-rate figures use every timepoint.']);
            rateWin = [-Inf, Inf];
        end
    elseif ischar(rateWindow)   % 'all'
        rateWin = [-Inf, Inf];
    else
        rateWin = sort(double(rateWindow(:)'));
    end
    fprintf('\nAreal rate vs dilatational flow, aligned times %g to %g min:\n', rateWin);

    % Keep each specimen's vertex values inside the window
    for gg = 1:numel(allRates)
        for ii = 1:numel(allRates{gg})
            r = allRates{gg}{ii};
            if isempty(r)
                continue
            end
            xAll{end + 1} = r.x; %#ok<AGROW>
            inTp = r.tau >= rateWin(1) & r.tau <= rateWin(2);
            inV  = inTp(r.it);
            r.x = r.x(inV);
            r.y = r.y(inV);
            r.timepoints = r.timepoints(inTp);
            r.nHoopsOut  = sum(r.nOutPerTp(inTp));
            r = rmfield(r, {'it', 'tau', 'tauRange', 'nOutPerTp'});
            if isempty(r.x)
                fprintf(2, '  %s: no vertex values inside the window, left out\n', r.name);
                failed(end+1, :) = {r.name, 'areal rate vs flow: no values inside rateWindow'}; %#ok<AGROW>
                r = [];
            else
                xWin{end + 1} = r.x; %#ok<AGROW>
                fprintf('  %s: %d timepoints, %d vertex values, %d outlier hoop(s) left out\n', ...
                    r.name, numel(r.timepoints), numel(r.x), r.nHoopsOut);
            end
            allRates{gg}{ii} = r;
        end
    end
end
clear rates

if all(cellfun(@isempty, [allRates{:}]))
    fprintf(2, '\nNo metric kinematics to use, so no areal rate vs dilatational flow figures.\n');
else
    % The same bins for every specimen and genotype, spanning the central 99%
    % of div(v_par). That range depends on which timepoints are pooled, so
    % both are printed: 'every timepoint' is what the earlier script used.
    rangeAll = centralRange(vertcat(xAll{:}));
    rangeWin = centralRange(vertcat(xWin{:}));
    clear xAll xWin
    fprintf(['\ndiv(v_par), central 99%%: [%.3g %.3g] 1/min at every timepoint on disk, ' ...
        '[%.3g %.3g] inside rateWindow\n'], rangeAll, rangeWin);
    if isempty(rateBinEdges)
        if strcmp(rateBinRange, 'all')
            rateBinEdges = linspace(rangeAll(1), rangeAll(2), nRateBins + 1);
            fprintf('Bins span the range at every timepoint (rateBinRange = ''all'').\n');
        else
            rateBinEdges = linspace(rangeWin(1), rangeWin(2), nRateBins + 1);
            fprintf('Bins span the range inside rateWindow (rateBinRange = ''window'').\n');
        end
    end
    rateCenters = (rateBinEdges(1:end-1) + rateBinEdges(2:end)) / 2;

    % Bin each specimen (its vertex values are dropped afterwards), then
    % average the specimens of each genotype
    rateAvg = cell(numel(genotypes), 1);
    yMax = 0;
    for gg = 1:numel(allRates)
        for ii = 1:numel(allRates{gg})
            if ~isempty(allRates{gg}{ii})
                r = allRates{gg}{ii};
                r.stats = binnedStats(r.x, r.y, rateBinEdges, minPointsPerBin);
                allRates{gg}{ii} = rmfield(r, {'x', 'y'});
                yMax = max(yMax, bandExtent(r.stats));
            end
        end
        have = ~cellfun(@isempty, allRates{gg});
        if any(have)
            rs = [allRates{gg}{have}];
            rateAvg{gg} = averageStats({rs.stats}, minSpecimensPerBin);
            rateAvg{gg}.names = {rs.name};
            yMax = max(yMax, bandExtent(rateAvg{gg}));
        end
    end

    % Axes shared by both figures: y symmetric about zero and, as in the
    % earlier figures, wide enough for every specimen's band too
    if isempty(rateXLim)
        rateXLim = rateBinEdges([1, end]);
    end
    if isempty(rateYLim)
        if yMax > 0
            rateYLim = 1.1 * yMax * [-1, 1];
        else
            rateYLim = [-0.1, 0.1];
        end
    end
    fprintf('\nAreal rate vs dilatational flow: %d bins, axes x = [%.3g %.3g], y = [%.3g %.3g] 1/min\n', ...
        numel(rateCenters), rateXLim, rateYLim);

    for gg = 1:numel(allRates)
        A = rateAvg{gg};
        if isempty(A)
            continue
        end
        if ~any(isfinite(A.mean))
            fprintf(2, 'No bin has %d %s specimens, so no areal-rate figure.\n', ...
                minSpecimensPerBin, genotypes(gg).name);
            continue
        end
        base = fullfile(outdir, ['arealRate_vs_divv_' genotypes(gg).name '_average']);
        ttl = sprintf('%s average, n = %d', genotypes(gg).title, numel(A.names));
        for ss = 1:numel(styles)
            fig = plotArealRate(rateCenters, A, styles{ss}, rateXLim, rateYLim, ttl, ...
                {'mean of specimens $\pm$ SEM', '$\pm$ 1 SD across specimens'});
            savePdfPng(fig, [base '_' styles{ss}]);
            close(fig);
        end
        fprintf('Saved %s_<style>.pdf/.png\n', base);
    end
end

%% -------------------------------------------------------------------------
%  Surface area and volume: each genotype's time-aligned ensemble
% --------------------------------------------------------------------------
ensGeno    = [];   % per ensemble: genotype index,
ensMembers = {};   % its specimens,
ensMinSpec = [];   % the fewest specimens a mean may use,
ensExtent  = [];   % and how far it has at least that many [min]
for gg = 1:numel(genotypes)
    gName   = genotypes(gg).name;
    members = allCurves{gg}(~cellfun(@isempty, allCurves{gg}));
    fprintf('\n%s ensemble, %d specimens:\n', gName, numel(members));
    for mm = 1:numel(members)
        c = members{mm};
        fprintf('  %s: lag %+g timepoints, aligned t0 = %d, data from %g to %g min\n', ...
            c.name, c.lag, c.t0, c.t(1), c.t(end));
    end
    if numel(members) < 2
        fprintf(2, '  fewer than 2 specimens to average, so no %s ensemble\n', gName);
        continue
    end
    if isempty(minSpecimensForMean)
        minSpec = numel(members);
    else
        minSpec = min(minSpecimensForMean, numel(members));
    end

    % How far this genotype has at least minSpec specimens
    E = alignEnsemble(members, cellfun(@(c) c.t0, members), minSpec, minutesPerTimepoint);
    fprintf('  all %d specimens have data from 0 to %g min; at least %d to %g min\n', ...
        numel(members), E.tau(find(E.nSpecimens == numel(members), 1, 'last')), ...
        minSpec, E.tau(end));
    ensGeno(end + 1)    = gg; %#ok<AGROW>
    ensMembers{end + 1} = members; %#ok<AGROW>
    ensMinSpec(end + 1) = minSpec; %#ok<AGROW>
    ensExtent(end + 1)  = E.tau(end); %#ok<AGROW>
end

% Both genotypes run from t = 0 to the same end time: the last time every
% genotype still has its minimum number of specimens (or ensembleEndTime,
% if earlier). At each time the mean and band use the specimens that have
% data there.
ens = {};   % one ensemble per genotype
endTime = [];
if ~isempty(ensMembers)
    [endTime, first] = min(ensExtent);
    reach = arrayfun(@(k) sprintf('%s has %d or more specimens to %g min', ...
        genotypes(ensGeno(k)).name, ensMinSpec(k), ensExtent(k)), 1:numel(ensExtent), ...
        'UniformOutput', false);
    if ~isempty(ensembleEndTime) && ensembleEndTime < endTime
        endTime = floor(ensembleEndTime / minutesPerTimepoint + 1e-9) * minutesPerTimepoint;
        fprintf('\nBoth panels run from t = 0 to %g min (ensembleEndTime).\n', endTime);
    else
        fprintf(['\nBoth panels run from t = 0 to %g min, the last time %s has at least ' ...
            '%d specimens.\n'], endTime, genotypes(ensGeno(first)).name, ensMinSpec(first));
    end
    fprintf('  (%s)\n', strjoin(reach, '; '));
    lastFrame = round(endTime / minutesPerTimepoint);
    for ee = 1:numel(ensMembers)
        gg      = ensGeno(ee);
        members = ensMembers{ee};
        E = alignEnsemble(members, cellfun(@(c) c.t0, members), ensMinSpec(ee), ...
            minutesPerTimepoint, lastFrame);
        E.genotype = genotypes(gg).name;
        E.title    = genotypes(gg).title;
        ens{end + 1} = E; %#ok<AGROW>
        fprintf('  %s: %s\n', E.genotype, coverageText(E));
    end
end

% Figures, with the same axes in both
if ~isempty(ens)
    if isempty(ensembleYLim)
        vals = [];
        for ee = 1:numel(ens)
            E  = ens{ee};
            bA = E.(['area' ensembleBand]);
            bV = E.(['volume' ensembleBand]);
            vals = [vals, E.areaMean + bA, E.areaMean - bA, ...
                E.volumeMean + bV, E.volumeMean - bV]; %#ok<AGROW>
        end
        ensembleYLim = paddedRange(vals(:));
    end
    ensembleXLim = [0, max(endTime, minutesPerTimepoint)];
    fprintf('\nSurface area and volume axes: x = [%g %g] min, y = [%.3g %.3g]\n', ...
        ensembleXLim, ensembleYLim);
    for ee = 1:numel(ens)
        E    = ens{ee};
        ttl  = sprintf('%s, n = %s, mean \\pm %s', E.title, nRangeText(E), ensembleBand);
        base = fullfile(outdir, ['ensemble_' E.genotype '_area_volume']);
        for ss = 1:numel(styles)
            fig = plotEnsemble(E, styles{ss}, ensembleXLim, ensembleYLim, ttl, ensembleBand);
            savePdfPng(fig, [base '_' styles{ss}]);
            close(fig);
        end
        fprintf('Saved %s_<style>.pdf/.png\n', base);
    end
end

%% -------------------------------------------------------------------------
%  Summary
% --------------------------------------------------------------------------
fprintf('\n========== Done ==========\n');
if isempty(failed)
    fprintf('All specimens processed. Results are in %s\n', outdir);
else
    fprintf(2, '%d problem(s):\n', size(failed, 1));
    for kk = 1:size(failed, 1)
        fprintf(2, '  %s: %s\n', failed{kk, 1}, failed{kk, 2});
    end
end

%% -------------------------------------------------------------------------
%  Local functions
% --------------------------------------------------------------------------

function dstr = specimenDate(dataDir)
% specimenDate  Acquisition date as 'YYYYMMDD', taken from the first date in
% the dataset path. Recognizes '20240528_...' and '2024-06-26_...' styles.
tok = regexp(dataDir, '(?<!\d)(20\d{2})-?(\d{2})-?(\d{2})(?!\d)', 'tokens', 'once');
if isempty(tok)
    error('No date (YYYYMMDD or YYYY-MM-DD) found in path:\n  %s', dataDir);
end
dstr = [tok{:}];
end

function TL = alignedTimeline(S, shiftSign)
% alignedTimeline  The shared aligned timeline of every dataset in the
% shifts file. Each dataset's lag onto the WT
% reference's timeline is its relative_shifts entry (empty or NaN for the
% reference, i.e. 0), plus cross_genotype_shift for bynMyo1c. t = 0 is the
% first aligned timepoint at which every dataset has data: the latest
% aligned start, from each dataset's first timepoint in its xp.mat.
canon = @(p) regexprep(regexprep(strtrim(char(p)), '[/\\]+', '/'), '/$', '');
n = numel(S.dataDirs);
TL.dataDirs  = cellfun(canon, S.dataDirs, 'UniformOutput', false);
TL.genotypes = S.genotypes;
TL.lags      = zeros(1, n);
TL.firsts    = nan(1, n);   % each dataset's first timepoint
for j = 1:n
    g = S.genotypes{j};
    p = find(S.([g '_indices']) == j);
    if isempty(p)
        error('Dataset %d (%s) is missing from %s_indices in the shifts file.', j, g, g);
    end
    sh = S.relative_shifts.(g);
    if iscell(sh)
        s = sh{p};
    else
        s = sh(p);
    end
    if ~isempty(s) && ~isnan(s)
        TL.lags(j) = double(s);
    end
    if ~strcmpi(g, 'wt')
        TL.lags(j) = TL.lags(j) + double(S.cross_genotype_shift);
    end
    try
        X = load(fullfile(S.dataDirs{j}, 'xp.mat'), 'xp');
        TL.firsts(j) = double(X.xp.fileMeta.timePoints(1));
    catch ME
        error('Cannot read the first timepoint of dataset %d (%s) from its xp.mat: %s', ...
            j, S.dataDirs{j}, ME.message);
    end
end
[TL.alignedTimepointAtZero, TL.lastToStart] = max(TL.firsts + shiftSign * TL.lags);
TL.shiftSign = shiftSign;
end

function k = findDataset(TL, dataDir)
% findDataset  Index of a dataset in the shifts file, found by its folder.
canon = @(p) regexprep(regexprep(strtrim(char(p)), '[/\\]+', '/'), '/$', '');
k = find(strcmp(TL.dataDirs, canon(dataDir)));
if isempty(k)
    error('%s is not in the shifts file''s dataDirs.', dataDir);
elseif numel(k) > 1
    error('%s is listed %d times in the shifts file''s dataDirs.', dataDir, numel(k));
end
end

function t0 = datasetT0(TL, k, shiftSign)
% datasetT0  Dataset k's aligned t0 on the timeline built with shiftSign:
% its timepoint that lands at t = 0, the first aligned timepoint at which
% every dataset in the shifts file has data,
%   t0 = alignedTimepointAtZero - shiftSign * lag,
% rounded to a whole timepoint when the lags are not
% whole numbers. Whether the dataset has data at t0 is up to the caller.
atZero = max(TL.firsts + shiftSign * TL.lags);
t0 = round(atZero - shiftSign * TL.lags(k));
end

function trimToAlignedMeshes(tubi)
% trimToAlignedMeshes  End tubi's timepoint list at its last aligned mesh.
%   Checks which timepoints listed in xp.mat have an aligned mesh on disk
%   and drops the trailing ones that don't. Errors if meshes are missing
%   before the last one. Only the in-memory TubULAR instance is changed,
%   never xp.mat.
tps  = tubi.xp.fileMeta.timePoints;
have = arrayfun(@(t) exist(sprintf(tubi.fullFileBase.alignedMesh, t), 'file') == 2, tps);
lastIdx = find(have, 1, 'last');
if isempty(lastIdx)
    error('No aligned meshes found in %s', tubi.dir.alignedMesh);
end

kept = tps(1:lastIdx);
gaps = kept(~have(1:lastIdx));
if ~isempty(gaps)
    error('Missing aligned meshes for t = %s(the last mesh on disk is t = %d)', ...
        sprintf('%d ', gaps), kept(end));
end

if lastIdx < numel(tps)
    fprintf('  Trimming timepoints to %d..%d (last aligned mesh); xp.mat lists up to %d\n', ...
        kept(1), kept(end), tps(end));
    tubi.xp.fileMeta.timePoints = kept;
    tubi.xp.tIdx = containers.Map(kept, 1:length(kept));   % keep TubULAR's index lookup in sync
end
end

function c = loadCurve(tubi, specOutdir, name, dataDir, lag, t0, mpt)
% loadCurve  One specimen's surface area and volume on the aligned
%   timeline: the raw values from measureSurfaceAreaVolume's cache, its
%   timepoints, and their aligned time (timepoint - t0) * mpt [min].
S   = load(fullfile(specOutdir, 'surfacearea_volume_stab.mat'), 'aas', 'vvs');
tps = tubi.xp.fileMeta.timePoints(:);
if numel(S.aas) ~= numel(tps)
    error(['Saved results have %d timepoints but the dataset has %d. ' ...
        'Rerun with overwrite = true.'], numel(S.aas), numel(tps));
end
if ~any(tps == t0)
    error('The aligned t0 = %g is not one of the timepoints, so the curves cannot be normalized.', t0);
end
c.name    = name;
c.dataDir = dataDir;
c.lag     = lag;                % timepoints onto the WT reference's timeline
c.t0      = t0;                 % aligned t0
c.tp      = tps;
c.t       = (tps - t0) * mpt;   % aligned time [min]
c.aRaw    = S.aas(:);
c.vRaw    = S.vvs(:);
end

function f = minutesPerUnit(units)
% minutesPerUnit  Factor that converts TubULAR's timeUnits to minutes.
switch lower(strtrim(units))
    case {'min', 'mins', 'minute', 'minutes'}
        f = 1;
    case {'s', 'sec', 'secs', 'second', 'seconds'}
        f = 1/60;
    case {'h', 'hr', 'hrs', 'hour', 'hours'}
        f = 60;
    otherwise
        warning('Unrecognized timeUnits ''%s''; treating them as minutes.', units);
        f = 1;
end
end

function lims = paddedRange(vals)
% paddedRange  [min max] of vals (NaNs ignored), widened by 5% on each side.
lo  = min(vals);
hi  = max(vals);
pad = 0.05 * (hi - lo);
if pad == 0
    pad = 0.5;
end
lims = [lo - pad, hi + pad];
end

function savePdfPng(fig, base)
% savePdfPng  Vector PDF plus 300 dpi PNG, with the page set to the figure
% size so all three versions come out the same size, panels in the same spot.
set(fig, 'Units', 'inches');
pos = get(fig, 'Position');
set(fig, 'PaperUnits', 'inches', 'PaperSize', pos(3:4), 'PaperPosition', [0 0 pos(3:4)]);
print(fig, [base '.pdf'], '-dpdf', '-painters');
print(fig, [base '.png'], '-dpng', '-r300');
end

function r = loadRates(tubi, mkRunFolder, name, specOutdir, t0, R)
% loadRates  One specimen's dilatational flow div(v_par) and areal rate of
% change 1/2 Tr[g^-1 gdot] (TubULAR's gdot, = div(v_par) - 2Hv_n) at every
% vertex in the middle of the AP axis, at every timepoint that has both
% files in the measurements folder measureMetricKinematics wrote. Hoops
% whose areal rate is far from typical are left out (the outlier check of
% the div(v_par) vs 2Hv_n script). Rates are converted to 1/min as in the
% strain-rate script. Each value keeps its timepoint: r.it indexes
% r.timepoints and r.tau, their aligned times (timepoint - t0) *
% minutesPerTimepoint, so the values can be cut to an aligned window.
if isempty(mkRunFolder)
    measDir = tubi.dir.metricKinematics.measurements;
else
    measDir = fullfile(tubi.dir.metricKinematics.root, mkRunFolder, 'measurements');
end
if ~exist(measDir, 'dir')
    error('Measurements folder not found: %s', measDir);
end
files = dir(fullfile(measDir, 'divv_vertices_*.mat'));
tps = [];
for kf = 1:numel(files)
    tp = sscanf(files(kf).name, 'divv_vertices_%d.mat');
    if ~isempty(tp) && exist(fullfile(measDir, sprintf('gdot_vertices_%06d.mat', tp)), 'file')
        tps(end + 1) = tp; %#ok<AGROW>
    end
end
if isempty(tps)
    error('No divv/gdot files in %s', measDir);
end
tps = sort(tps);
nT = numel(tps);

% [AP x around x time] over the middle of the AP axis, without the
% duplicate seam column
for kk = 1:nT
    tmp = load(fullfile(measDir, sprintf('divv_vertices_%06d.mat', tps(kk))), R.dVar);
    dd  = tmp.(R.dVar);
    tmp = load(fullfile(measDir, sprintf('gdot_vertices_%06d.mat', tps(kk))), R.gVar);
    gd  = tmp.(R.gVar);
    if kk == 1
        [nU, nV] = size(dd);
        padN   = round(R.apTrim * nU);
        apKeep = (padN + 1):(nU - padN);
        D = nan(numel(apKeep), nV - 1, nT, 'single');
        G = D;
    end
    D(:, :, kk) = dd(apKeep, 1:nV-1);
    G(:, :, kk) = gd(apKeep, 1:nV-1);
end
rateScale = tubi.timeInterval / R.minutesPerTimepoint;   % saved rates -> 1/min
D = D * rateScale;
G = G * rateScale;

% Leave out hoops (one AP position at one timepoint) whose areal rate is
% far from typical somewhere around the tube
gv = double(G(isfinite(G)));
gCenter = median(gv);
gSD = 1.4826 * median(abs(gv - gCenter));   % robust SD (MAD scaled to match an SD)
if ~(gSD > 0)
    gSD = eps;
end
dev = abs(G - gCenter);
hoopOut = any(dev > R.zThresh * gSD & dev > R.minDeviation, 2);   % [AP x 1 x time]
keep = ~repmat(hoopOut, [1, size(G, 2), 1]) & isfinite(D) & isfinite(G);
it = repmat(reshape(uint32(1:nT), 1, 1, nT), size(G, 1), size(G, 2));   % timepoint index

r.name       = name;
r.specOutdir = specOutdir;
r.measDir    = measDir;
r.t0         = t0;
r.timepoints = tps;
r.tau        = (tps - t0) * R.minutesPerTimepoint;   % aligned time of each timepoint [min]
r.tauRange   = [min(r.tau), max(r.tau)];
r.nOutPerTp  = reshape(sum(hoopOut, 1), 1, []);      % outlier hoops at each timepoint
r.x  = D(keep);
r.y  = G(keep);
r.it = it(keep);
fprintf('  areal rate vs flow: %d timepoints on disk, aligned t = %g to %g min\n', ...
    nT, r.tauRange);
end

function lims = centralRange(v)
% centralRange  [lo hi] of the central 99% of v: its values at ranks 0.5%
% and 99.5%, which is how the areal-rate bins are set (as in the earlier
% script).
v = sort(double(v(:)));
n = numel(v);
if n == 0
    lims = [NaN, NaN];
    return
end
lims = [v(max(1, round(0.005 * n))), v(max(1, round(0.995 * n)))];
end

function s = binnedStats(x, y, edges, minN)
% binnedStats  Count, mean, SD and SEM of y within bins of x; bins with
% fewer than minN points are NaN.
nb = numel(edges) - 1;
b  = discretize(double(x(:)), edges);
ok = ~isnan(b);
b  = b(ok);
yy = double(y(:));
yy = yy(ok);
s.n    = accumarray(b, 1, [nb, 1])';
s.mean = accumarray(b, yy, [nb, 1], @mean, NaN)';
s.sd   = accumarray(b, yy, [nb, 1], @std, NaN)';
s.sem  = s.sd ./ sqrt(s.n);
bad = s.n < minN;
s.mean(bad) = NaN;
s.sd(bad)   = NaN;
s.sem(bad)  = NaN;
end

function s = averageStats(statsList, minSpec)
% averageStats  A genotype's average: per bin, the mean of its specimens'
% bin means, with SD and SEM across specimens. Bins with fewer than minSpec
% specimens are NaN.
M = cell2mat(cellfun(@(st) st.mean(:)', statsList(:), 'UniformOutput', false));
s.nSpec = sum(isfinite(M), 1);
s.mean  = mean(M, 1, 'omitnan');
s.sd    = std(M, 0, 1, 'omitnan');
s.sem   = s.sd ./ sqrt(s.nSpec);
bad = s.nSpec < minSpec;
s.mean(bad) = NaN;
s.sd(bad)   = NaN;
s.sem(bad)  = NaN;
s.specimenMeans = M;
end

function m = bandExtent(st)
% bandExtent  Largest |mean +/- SD| over the bins that have data (0 if none).
v = abs([st.mean + st.sd, st.mean - st.sd]);
v = v(isfinite(v));
if isempty(v)
    m = 0;
else
    m = max(v);
end
end

function fig = plotArealRate(centers, st, style, xLims, yLims, titleStr, legendNames)
% plotArealRate  Binned areal rate of change against dilatational flow, as
% in TubULAR's Extended Data Fig. 4d: the mean per bin (black dots and
% line, error bars = SEM) over a gray band of +/- 1 SD, with dashed lines
% at zero. style is 'labeled', 'ticks' or 'bare'; the axes sit at the same
% place in all three.

% Layout, in inches
panelW = 2.6;   panelH = 2.0;
padL   = 1.2;   padR   = 0.25;
padT   = 0.45;  padB   = 0.7;
figW = padL + panelW + padR;
figH = padB + panelH + padT;
fig  = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [1 1 figW figH]);
ax   = axes('Parent', fig, 'Units', 'inches', 'Position', [padL, padB, panelW, panelH]);
hold(ax, 'on')

zeroColor = [0.3, 0.6, 0.9];
ok = isfinite(st.mean);
c  = centers(ok);
m  = st.mean(ok);
sd = st.sd(ok);
if any(ok)
    hBand = fill(ax, [c, fliplr(c)], [m + sd, fliplr(m - sd)], [0.8, 0.8, 0.8], 'EdgeColor', 'none');
end
plot(ax, [0, 0], yLims, '--', 'Color', zeroColor, 'LineWidth', 0.75)
plot(ax, xLims, [0, 0], '--', 'Color', zeroColor, 'LineWidth', 0.75)
if any(ok)
    hMean = errorbar(ax, c, m, st.sem(ok), 'k.-', 'MarkerSize', 10, 'LineWidth', 1, 'CapSize', 4);
end
set(ax, 'XLim', xLims, 'YLim', yLims, 'Box', 'on', 'TickDir', 'in', 'FontSize', 8, 'LineWidth', 0.5)
if strcmp(style, 'bare')
    set(ax, 'XTickLabel', [], 'YTickLabel', [])
end
if strcmp(style, 'labeled')
    xlabel(ax, 'dilatational flow, $\nabla\cdot\mathbf{v}_\parallel$ [1/min]', 'Interpreter', 'latex')
    ylabel(ax, {'areal rate of change,', ...
        '$\frac{1}{2}\mathrm{Tr}\left[\mathbf{g}^{-1}\dot{\mathbf{g}}\right]$ [1/min]'}, ...
        'Interpreter', 'latex')
    title(ax, titleStr, 'Interpreter', 'tex', 'FontWeight', 'normal', 'FontSize', 8)
    if any(ok)
        legend(ax, [hMean, hBand], legendNames, 'Interpreter', 'latex', 'Location', 'northwest', ...
            'Box', 'off', 'FontSize', 7)
    end
end
end

function E = alignEnsemble(members, t0s, minSpec, mpt, lastFrame)
% alignEnsemble  Specimens on the aligned timeline: each one's frames are
% counted from its aligned t0 (frame = timepoint - t0), so frame 0 is t = 0
% for all of them, and its area and volume are normalized to their value
% there. The ensemble runs from t = 0 to frame lastFrame (default: the last
% frame at which at least minSpec specimens have data). At each time the
% mean, SD and SEM are over the specimens that have data there, and NaN
% where fewer than minSpec do.
n = numel(members);
minSpec = min(minSpec, n);
names = cellfun(@(c) c.name, members, 'UniformOutput', false);
fr = cell(1, n);   % frames relative to each specimen's t0
i0 = zeros(1, n);  % index of each specimen's t0
for mm = 1:n
    fr{mm} = members{mm}.tp(:)' - t0s(mm);
    hit = find(fr{mm} == 0, 1);
    if ~isempty(hit)
        i0(mm) = hit;
    end
end
if any(i0 == 0)
    error('t = 0 is not one of the timepoints of %s.', strjoin(names(i0 == 0), ', '));
end
if nargin < 5 || isempty(lastFrame)
    starts = cellfun(@min, fr);
    ends   = cellfun(@max, fr);
    frames = 0:max(ends);
    nCover = sum(frames >= starts(:) & frames <= ends(:), 1);   % specimens covering each time
    lastFrame = frames(find(nCover >= minSpec, 1, 'last'));
end
frames = 0:lastFrame;
A = nan(n, numel(frames));
V = nan(n, numel(frames));
for mm = 1:n
    c = members{mm};
    a = c.aRaw(:)';
    v = c.vRaw(:)';
    [inRange, loc] = ismember(fr{mm}, frames);
    A(mm, loc(inRange)) = a(inRange) / a(i0(mm));
    V(mm, loc(inRange)) = v(inRange) / v(i0(mm));
end
nAV = sum(isfinite(A) & isfinite(V), 1);   % specimens with area and volume at each time
E.names      = names;
E.dataDirs   = cellfun(@(c) c.dataDir, members, 'UniformOutput', false);
E.lags       = cellfun(@(c) c.lag, members);
E.t0         = t0s(:)';
E.minutesPerTimepoint = mpt;
E.minSpecimens = minSpec;
E.tau        = frames * mpt;
E.area       = A;
E.volume     = V;
E.nSpecimens = nAV;
E.areaMean   = mean(A, 1, 'omitnan');
E.areaSD     = std(A, 0, 1, 'omitnan');
E.areaSEM    = E.areaSD ./ sqrt(nAV);
E.volumeMean = mean(V, 1, 'omitnan');
E.volumeSD   = std(V, 0, 1, 'omitnan');
E.volumeSEM  = E.volumeSD ./ sqrt(nAV);

% No mean or band where too few specimens have data
few = nAV < minSpec;
for q = {'area', 'volume'}
    E.([q{1} 'Mean'])(few) = NaN;
    E.([q{1} 'SD'])(few)   = NaN;
    E.([q{1} 'SEM'])(few)  = NaN;
end
end

function txt = coverageText(E)
% coverageText  How many specimens the ensemble mean uses over time, e.g.
% 'n = 5 to 84 min, 4 to 98 min, 3 to 112 min'.
n = E.nSpecimens;
runEnd = [find(diff(n) ~= 0), numel(n)];   % last frame of each stretch with the same n
parts = {};
for k = 1:numel(runEnd)
    nk = n(runEnd(k));
    if nk < E.minSpecimens
        parts{end + 1} = sprintf('then no mean (fewer than %d specimens) to %g min', ...
            E.minSpecimens, E.tau(end)); %#ok<AGROW>
        break
    end
    parts{end + 1} = sprintf('%d to %g min', nk, E.tau(runEnd(k))); %#ok<AGROW>
end
txt = ['n = ' strjoin(parts, ', ')];
end

function s = nRangeText(E)
% nRangeText  The number of specimens behind the plotted mean, e.g. '5' or
% '3-5' when specimens drop out before the end.
n = E.nSpecimens(E.nSpecimens >= E.minSpecimens);
if isempty(n)
    s = '0';
elseif min(n) == max(n)
    s = sprintf('%d', max(n));
else
    s = sprintf('%d-%d', min(n), max(n));
end
end

function fig = plotEnsemble(E, style, xLims, yLims, titleStr, band)
% plotEnsemble  Time-aligned ensemble mean of the normalized surface area
% and volume, each with a shaded band (SD or SEM across specimens). style
% is 'labeled', 'ticks' or 'bare'; the axes sit at the same place in all
% three.

% Layout, in inches
panelW = 3.2;   panelH = 2.2;
padL   = 0.9;   padR   = 0.25;
padT   = 0.45;  padB   = 0.65;
figW = padL + panelW + padR;
figH = padB + panelH + padT;
fig  = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [1 1 figW figH]);
ax   = axes('Parent', fig, 'Units', 'inches', 'Position', [padL, padB, panelW, panelH]);
hold(ax, 'on')

q      = {'area', 'volume'};
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980]};   % surface area blue, volume orange
hLine  = gobjects(1, numel(q));
for kk = 1:numel(q)
    m  = E.([q{kk} 'Mean']);
    b  = E.([q{kk} band]);
    ok = isfinite(m) & isfinite(b);
    t  = E.tau(ok);
    fill(ax, [t, fliplr(t)], [m(ok) + b(ok), fliplr(m(ok) - b(ok))], colors{kk}, ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none');
    hLine(kk) = plot(ax, E.tau, m, '-', 'Color', colors{kk}, 'LineWidth', 1.5);
end
set(ax, 'XLim', xLims, 'YLim', yLims, 'Box', 'on', 'TickDir', 'in', 'FontSize', 8, 'LineWidth', 0.5)
if strcmp(style, 'bare')
    set(ax, 'XTickLabel', [], 'YTickLabel', [])
end
if strcmp(style, 'labeled')
    xlabel(ax, 'aligned time [min]')
    ylabel(ax, 'normalized area or volume')
    title(ax, titleStr, 'Interpreter', 'tex', 'FontWeight', 'normal', 'FontSize', 8)
    legend(ax, hLine, {'surface area', 'volume'}, 'Location', 'best', 'Box', 'off', 'FontSize', 7)
end
end
