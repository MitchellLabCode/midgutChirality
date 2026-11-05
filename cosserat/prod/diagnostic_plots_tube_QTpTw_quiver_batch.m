% Batch diagnostic tube + material-frame plots for the (T, tw, p) grid.
%
% Same figure as diagnostic_plots_tube_QTpTw_batch.m -- tube renderings of the
% reference and simulated centerlines over a sampled (T, tw, p) grid -- with the
% material frame overlaid as quivers, in the style of plot_tube_quiver.m /
% plot_3D_quiver.m. Reference directors are drawn in black, simulated ones in
% the sample's vik colour, so a twist mismatch is visible directly rather than
% only through the centerline shape.
%
% Sampling: you list a few values for each of T, tw and p; each is snapped to the
% nearest value actually present on disk. All sampled T are shown on the same
% image: rows = (T, tw) blocks, cols = p (one figure per sample, per view).
%
% Colour convention (simulated tube + its quivers), matching the other batch
% scripts:
%   even-numbered samples (s2, s4, ...) -> WT  colour  vik(22,:)
%   odd-numbered  samples (s1, s3, ...) -> Myo colour  vik(234,:)
%
% Two deliberate differences from plot_tube_quiver.m, both because this is a
% multi-panel comparison rather than a single large frame:
%
%   * Arrows are drawn at a FIXED length in microns (AutoScale off).  With
%     'AutoScale','on' MATLAB normalizes arrow length per axes, so the same
%     director would render at different lengths in different tiles and the
%     panels could not be compared.  See quiver_len below.
%   * Only every quiver_stride-th director is drawn.  The reference frame has
%     ~99 elements and the simulation ~49; all of them in a 220x200 px tile is
%     an unreadable thicket.
%
% No local functions are used, so this can be pasted into the command window or
% run on older MATLAB releases.

addpath(".")

% ======================= CONFIG ============================================
% Folders to compare. Leave the date as * to pick up every QTpTw run, or pin a
% date_str (e.g. "20260723") to a single sweep.
folder_glob = char("*_QTpTw_contact2_t*_w3_a3e2_s*V6");
%folder_glob = char("*_QTpTw_contact2_t*_w3_a3e2_s7V6");

% Representative grid points to show (snapped to nearest available on disk).
% NB: the plastic grid of the FINE sweep is 10.^(-2:0.3:2), whose exponent steps
% of 0.3 never land on 0 or 2 -- so there p = 1 and p = 100 are off-grid and snap
% to ~1.26 and ~79.4. The coarse (20260723) sweeps do contain 1 and 100. Any
% off-grid request is snapped and reported.
T_samp  = [ 3.0];      % simulation time
tw_samp = [ 1.0, ];      % twist-to-bend ratio B3/B1
p_samp  = [ 1.0, 10, ];   % plasticity (viscoelastic) rate

% To get the B_t/B = 1 sheet over the full (T, p) grid of a coarse sweep, use:
%   T_samp  = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0];
%   tw_samp = [1.0];
%   p_samp  = [0.0, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0];

views_to_plot = {'ventral', 'side'};   % any of {'ventral','side','dorsal'}
% 'twin'     -> every view is drawn side by side WITHIN each (T,tw,p) cell, so a
%               single figure shows both; matches the _twin_ convention of
%               diagnostic_twp_quiver.m.  Columns become p x view.
% 'separate' -> one figure per view (the old behaviour).
view_mode = 'twin';

tube_r = 5;                      % tube radius (microns)
tube_n = 24;                     % circumferential resolution
save_dpi = '-r500';

% ---- material-frame quivers ----
draw_quivers  = true;            % set false to fall back to tubes only
director_idx  = 1;               % which director to draw: 1 = d1, 2 = d2, 3 = d3
quiver_len    = 25;              % arrow length in MICRONS (fixed; AutoScale off)
quiver_stride = 6;               % draw every Nth director along the rod
quiver_lw     = 1.0;             % arrow line width
draw_ref_frame = true;           % also draw the reference (imaged) frame in black
% ===========================================================================

% ---- colours ---------------------------------------------------------------
vik = load('vik.mat').vik;
wt_color  = vik(22,:);   % WT
myo_color = vik(234,:);  % Myo

% ---- discover output folders -----------------------------------------------
listing = dir(folder_glob);
listing = listing([listing.isdir]);   % directories only
if isempty(listing)
    error('No output folders matching "%s" found in %s', folder_glob, pwd);
end

% ---- shared axes limits ----------------------------------------------------
% ONE bounding box across every folder processed in this run, so panels are
% comparable within a figure AND figures are comparable with each other.
% num_t is read with matfile() so the reference array's size is available
% without loading it, matching the t0:num_t range the plotting loop uses.
xmin =  1e9; xmax = -1e9;
ymin =  1e9; ymax = -1e9;
zmin =  1e9; zmax = -1e9;
for f = 1:numel(listing)
    tokb = regexp(listing(f).name, '_t(\d+)_.*_s(\d+)V6', 'tokens');
    if isempty(tokb)
        continue
    end
    t0b  = str2double(tokb{1}{1});
    sidb = str2double(tokb{1}{2});
    pmb  = sprintf('../data/meshes_spRS/meshes_sample%d.mat', sidb);
    pdb  = sprintf('../data/refIV6_sample%d_walls3_cball.mat', sidb);
    if ~isfile(pmb) || ~isfile(pdb)
        continue
    end
    ntb = size(matfile(pdb), 'x', 1);
    mb  = load(pmb); mb = mb.meshes;
    for t = t0b:min(ntb, numel(mb))
        v = mb{t}.v;
        xmin = min(xmin, min(v(:,1)));  xmax = max(xmax, max(v(:,1)));
        ymin = min(ymin, min(v(:,2)));  ymax = max(ymax, max(v(:,2)));
        zmin = min(zmin, min(v(:,3)));  zmax = max(zmax, max(v(:,3)));
    end
    clear mb
end
if xmin > xmax
    error('Could not build a bounding box from the folders matched by "%s".', folder_glob);
end
fprintf('shared axes limits: x [%.0f %.0f]  y [%.0f %.0f]  z [%.0f %.0f]\n', ...
        xmin, xmax, ymin, ymax, zmin, zmax);

for f = 1:numel(listing)
    folder_name = listing(f).name;

    % parse t0 and sample number from the folder name
    tok = regexp(folder_name, '_t(\d+)_.*_s(\d+)V6', 'tokens');
    if isempty(tok)
        fprintf('Skipping unrecognised folder: %s\n', folder_name);
        continue
    end
    t0        = str2double(tok{1}{1});
    sample_id = str2double(tok{1}{2});

    % colour by parity of the sample number
    if mod(sample_id, 2) == 0
        vik_c = wt_color;  color_tag = 'WT';
    else
        vik_c = myo_color; color_tag = 'Myo';
    end

    % matching mesh + reference data for this sample
    path_meshes = sprintf('../data/meshes_spRS/meshes_sample%d.mat', sample_id);
    path_data   = sprintf('../data/refIV6_sample%d_walls3_cball.mat', sample_id);
    if ~isfile(path_meshes) || ~isfile(path_data)
        fprintf('Skipping sample %d: missing mesh or data file.\n', sample_id);
        continue
    end

    fprintf('==== sample %d (%s), t0=%d, folder=%s ====\n', ...
            sample_id, color_tag, t0, folder_name);

    meshes = load(path_meshes); meshes = meshes.meshes;
    data   = load(path_data);
    data_x = data.x;
    data_Q = data.Q;               % reference material frames (num_t, N-1, 3, 3)
    num_t  = size(data_x, 1);

    % axes limits come from the shared bounding box computed above
    verts = meshes{num_t}.v;
    faces = meshes{num_t}.f;

    % rescale simulated centerline (non-dimensional) back to microns by L(t0)
    L_t0 = sum(sqrt(sum(diff(squeeze(data_x(t0,:,:)), 1, 1).^2, 2)));

    % reference frame at the final time, subsampled
    d1_ref  = squeeze(data_Q(num_t, :, director_idx, :));      % (N-1) x 3
    ref_sub = 1:quiver_stride:size(d1_ref, 1);
    ref_bx  = squeeze(data_x(num_t, ref_sub, 1));
    ref_by  = squeeze(data_x(num_t, ref_sub, 2));
    ref_bz  = squeeze(data_x(num_t, ref_sub, 3));

    % index the available (T, tw, p) files
    mat_files = dir(fullfile(folder_name, 'cutoff01contact_T*_tw*_p*.mat'));
    if isempty(mat_files)
        fprintf('Skipping sample %d: no sweep .mat files in %s\n', sample_id, folder_name);
        continue
    end
    fileT = nan(numel(mat_files),1);
    fileW = nan(numel(mat_files),1);
    fileP = nan(numel(mat_files),1);
    for k = 1:numel(mat_files)
        tk = regexp(mat_files(k).name, 'T([\d.]+)_tw([\d.]+)_p([\d.]+)\.mat', 'tokens');
        if ~isempty(tk)
            fileT(k) = str2double(tk{1}{1});
            fileW(k) = str2double(tk{1}{2});
            fileP(k) = str2double(tk{1}{3});
        end
    end
    keep = ~isnan(fileT);
    mat_files = mat_files(keep); fileT = fileT(keep); fileW = fileW(keep); fileP = fileP(keep);

    % snap requested values to nearest available
    snap = @(req, avail) arrayfun(@(vv) avail(find(abs(avail - vv) == min(abs(avail - vv)), 1)), req);
    T_use  = snap(T_samp,  unique(fileT));
    tw_use = snap(tw_samp, unique(fileW));
    p_use  = snap(p_samp,  unique(fileP));

    % report any value that had to move (inlined; no local function)
    reqs = {T_samp, tw_samp, p_samp};  uses = {T_use, tw_use, p_use};  nms = {'T','tw','p'};
    for a = 1:3
        for i = 1:numel(reqs{a})
            if abs(reqs{a}(i) - uses{a}(i)) > 1e-9
                fprintf('  %s: requested %g -> nearest available %g\n', nms{a}, reqs{a}(i), uses{a}(i));
            end
        end
    end

    reltol = 1e-6;
    n_rows = numel(T_use) * numel(tw_use);   % rows = (T, tw) blocks

    if strcmp(view_mode, 'twin')
        figure_specs = {views_to_plot};      % one figure, all views per cell
    else
        figure_specs = cell(1, numel(views_to_plot));
        for iv = 1:numel(views_to_plot)
            figure_specs{iv} = views_to_plot(iv);
        end
    end

    % all sampled T stacked on the same image
    for iv = 1:numel(figure_specs)
        cell_views = figure_specs{iv};       % views drawn inside each cell
        n_cols = numel(p_use) * numel(cell_views);   % cols = p x view
        view_tag = strjoin(cell_views, '_');

        fig = figure('Position', [10, 10, 220*n_cols+120, 200*n_rows+120]);
        tiledlayout(n_rows, n_cols, 'Padding', 'tight', 'TileSpacing', 'compact');

        for it = 1:numel(T_use)               % row-blocks = simulation time
            Tval = T_use(it);
            for iw = 1:numel(tw_use)          % rows within block = twist ratio
                Wval = tw_use(iw);
                for ip = 1:numel(p_use)       % cols = plasticity (x view)
                    Pval = p_use(ip);
                    % NB: no nexttile here -- tiles are claimed inside the
                    % cell_views loop below, one per view of this (T,tw,p) cell

                    match = find(abs(fileT - Tval) <= reltol*max(1,abs(Tval)) & ...
                                 abs(fileW - Wval) <= reltol*max(1,abs(Wval)) & ...
                                 abs(fileP - Pval) <= reltol*max(1,abs(Pval)), 1);

                    if isempty(match)
                        for icv = 1:numel(cell_views)
                            nexttile;
                            title(sprintf('T=%g, tw=%g, p=%g\n(missing)', Tval, Wval, Pval), 'fontsize', 6);
                            axis off
                        end
                        continue
                    end

                    sim_data = load(fullfile(folder_name, mat_files(match).name));
                    xsaves = sim_data.xsaves;
                    Qsaves = sim_data.Qsaves;

                    for icv = 1:numel(cell_views)   % one tile per view of this cell
                        this_view = cell_views{icv};
                        nexttile;

                        trisurf(triangulation(faces, verts), 'edgecolor', 'none', ...
                                'facecolor', 0.9*[1,1,1], 'faceAlpha', 0.2);
                        hold on

                        [xx, yy, zz] = tubeplot(squeeze(data_x(num_t,:,:))', tube_r, tube_n);
                        surf(xx, yy, zz, 'Facecolor', 0.6*[1,1,1], 'EdgeColor', 'none');

                        sim_x = [squeeze(xsaves(end,:,1)*L_t0); ...
                                 squeeze(xsaves(end,:,2)*L_t0); ...
                                 squeeze(xsaves(end,:,3)*L_t0)];
                        [xx, yy, zz] = tubeplot(sim_x, tube_r, tube_n);
                        surf(xx, yy, zz, 'Facecolor', vik_c, 'EdgeColor', 'none');

                        % ---- material-frame quivers -------------------------------
                        if draw_quivers
                            % simulated frame: bases are the element start nodes, so
                            % drop the last centerline node to match Qsaves
                            d1_sim  = squeeze(Qsaves(end, :, director_idx, :));   % Ns x 3
                            sim_sub = 1:quiver_stride:size(d1_sim, 1);
                            sbx = squeeze(xsaves(end, sim_sub, 1)) * L_t0;
                            sby = squeeze(xsaves(end, sim_sub, 2)) * L_t0;
                            sbz = squeeze(xsaves(end, sim_sub, 3)) * L_t0;

                            if draw_ref_frame
                                quiver3(ref_bx(:), ref_by(:), ref_bz(:), ...
                                        d1_ref(ref_sub,1)*quiver_len, ...
                                        d1_ref(ref_sub,2)*quiver_len, ...
                                        d1_ref(ref_sub,3)*quiver_len, ...
                                        0, 'color', 'k', 'linewidth', quiver_lw);
                            end
                            quiver3(sbx(:), sby(:), sbz(:), ...
                                    d1_sim(sim_sub,1)*quiver_len, ...
                                    d1_sim(sim_sub,2)*quiver_len, ...
                                    d1_sim(sim_sub,3)*quiver_len, ...
                                    0, 'color', vik_c, 'linewidth', quiver_lw);
                            % the trailing 0 disables MATLAB's per-axes autoscaling,
                            % so quiver_len is honoured identically in every tile
                        end

                        axis([xmin, xmax, ymin, ymax, zmin, zmax])
                        daspect([1,1,1])
                        axis off
                        % only the first tile of a cell repeats the parameters, so
                        % the twin layout does not print them twice
                        if icv == 1
                            title(sprintf('T=%g, tw=%g, p=%g\n%s', Tval, Wval, Pval, ...
                                          this_view), 'fontsize', 6);
                        else
                            title(this_view, 'fontsize', 6);
                        end

                        lighting gouraud
                        material dull
                        if strcmp(this_view, 'side')
                            view([0, 0]);    % side view
                        elseif strcmp(this_view, 'dorsal')
                            view(0, 90);     % dorsal view
                        else
                            view(0, -90);    % ventral view
                        end
                        camlight('headlight');
                    end   % cell_views
                end
            end
        end

        if draw_quivers
            qtag = sprintf('d%d frame: black = reference, colour = simulation', director_idx);
        else
            qtag = 'tubes only';
        end
        sgtitle(sprintf('sample %d (%s), t0=%d,  %s  (rows: T x tw, cols: p x view)\n%s', ...
                        sample_id, color_tag, t0, strjoin(cell_views, ' + '), qtag), ...
                'fontsize', 9);

        save_filename = sprintf('%s_tubegridq_allT_%s.png', char(folder_name), view_tag);
        print(save_filename, '-dpng', save_dpi);
        fprintf('Saved %s\n', save_filename);
        close(fig)
    end
end

fprintf('~done~\n');
