% Timelapse tube plots at ONE (T, tw, p) operating point, for one or more samples.
%
% Companion to diagnostic_plots_tube_QTpTw_quiver_batch.m, which sweeps
% parameters at the final time. This one fixes the parameters and sweeps TIME:
% rows = views (ventral / side), cols = time points, each panel showing the
% imaged rod (grey) and the simulation (sample colour) together with their
% material frames, so the two can be followed as they diverge.
%
% List several ids in samples_to_plot to get one figure per sample in a single
% run; each sample's output folder is found by glob, so only the ids are needed.
%
% Same rendering settings as the grid script: tube radius/resolution, fixed
% quiver length in microns, quiver stride, vik colour by sample parity.
%
% The simulation's n_saves snapshots span [0, T], which corresponds to
% reference frames t0..num_t, so the simulation is interpolated onto the
% requested reference frames the way plot_tube_quiver.m does.
%
% No local functions are used, so this can be pasted into the command window.

addpath(".")

% ======================= CONFIG ============================================
samples_to_plot = [7];           % e.g. [4 5 6] -- one figure per sample

% Each sample's output folder is found with this pattern (%d = sample id). If
% several match, the last by name (the most recent date stamp) is used.
folder_pattern = '*_QTpTw_contact2_t*_s%dV6';

% target operating point (snapped to the nearest file actually on disk)
T_target  = 3.0;
tw_target = 1.0;
p_target  = 10.0;

% time points, as REFERENCE FRAME indices. Leave empty to use n_times frames
% evenly spaced from t0 to the last frame.
times_to_show = [];
n_times       = 5;

dt_min = 2;                      % minutes between reference frames

views_to_plot = {'ventral', 'side'};   % rows; any of ventral / side / dorsal

tube_r = 5;                      % tube radius (microns)
tube_n = 24;                     % circumferential resolution
save_dpi = '-r400';

% ---- material-frame quivers (same convention as the grid script) ----
draw_quivers   = true;
director_idx   = 1;              % 1 = d1, 2 = d2, 3 = d3
quiver_len     = 25;             % arrow length in MICRONS (fixed; AutoScale off)
quiver_stride  = 6;              % draw every Nth director
quiver_lw      = 1.0;
draw_ref_frame = true;           % imaged frame in black alongside the simulated one
% ===========================================================================

vik = load('vik.mat').vik;

% ---- shared axes limits ----------------------------------------------------
% ONE bounding box for every panel of every sample, so tiles are comparable
% within a figure and figures are comparable with each other. Accumulated over
% each sample's meshes from its own t0 (parsed from the folder name) onwards.
xmin =  1e9; xmax = -1e9; ymin =  1e9; ymax = -1e9; zmin =  1e9; zmax = -1e9;
for isamp = 1:numel(samples_to_plot)
    sid = samples_to_plot(isamp);
    pm  = sprintf('../data/meshes_spRS/meshes_sample%d.mat', sid);
    cnd = dir(sprintf(folder_pattern, sid));  cnd = cnd([cnd.isdir]);
    if ~isfile(pm) || isempty(cnd)
        continue
    end
    [~, o] = sort({cnd.name});
    tk = regexp(cnd(o(end)).name, '_t(\d+)_', 'tokens');
    if isempty(tk)
        continue
    end
    t0b = str2double(tk{1}{1});
    mb  = load(pm); mb = mb.meshes;
    for t = t0b:numel(mb)
        v = mb{t}.v;
        xmin = min(xmin, min(v(:,1)));  xmax = max(xmax, max(v(:,1)));
        ymin = min(ymin, min(v(:,2)));  ymax = max(ymax, max(v(:,2)));
        zmin = min(zmin, min(v(:,3)));  zmax = max(zmax, max(v(:,3)));
    end
    clear mb
end
if xmin > xmax
    error('Could not build a bounding box: no usable sample in samples_to_plot.');
end
fprintf('shared axes limits: x [%.0f %.0f]  y [%.0f %.0f]  z [%.0f %.0f]\n', ...
        xmin, xmax, ymin, ymax, zmin, zmax);

for isamp = 1:numel(samples_to_plot)
    sample_id = samples_to_plot(isamp);

    % ---- locate this sample's output folder -------------------------------
    cand = dir(sprintf(folder_pattern, sample_id));
    cand = cand([cand.isdir]);
    if isempty(cand)
        fprintf('Skipping sample %d: no folder matching "%s"\n', ...
                sample_id, sprintf(folder_pattern, sample_id));
        continue
    end
    [~, ord] = sort({cand.name});
    folder_name = cand(ord(end)).name;      % most recent date stamp
    if numel(cand) > 1
        fprintf('sample %d: %d folders match, using %s\n', ...
                sample_id, numel(cand), folder_name);
    end

% ---- colours (even sample -> WT, odd -> Myo), matching the batch scripts ---
if mod(sample_id, 2) == 0
    vik_c = vik(22,:);   color_tag = 'WT';
else
    vik_c = vik(234,:);  color_tag = 'Myo';
end

% ---- reference data + meshes ----------------------------------------------
path_meshes = sprintf('../data/meshes_spRS/meshes_sample%d.mat', sample_id);
path_data   = sprintf('../data/refIV6_sample%d_walls3_cball.mat', sample_id);
if ~isfile(path_meshes) || ~isfile(path_data)
    fprintf('Skipping sample %d: missing mesh or reference data.\n', sample_id);
    continue
end
meshes = load(path_meshes); meshes = meshes.meshes;
data   = load(path_data);
data_x = data.x;
data_Q = data.Q;
num_t  = size(data_x, 1);

% t0 from the folder name (same convention as the batch scripts)
tok = regexp(char(folder_name), '_t(\d+)_', 'tokens');
if isempty(tok)
    fprintf('Skipping sample %d: cannot parse t0 from "%s".\n', sample_id, folder_name);
    continue
end
t0 = str2double(tok{1}{1});

% ---- locate the file for the requested (T, tw, p) -------------------------
mat_files = dir(fullfile(folder_name, 'cutoff01contact_T*_tw*_p*.mat'));
if isempty(mat_files)
    fprintf('Skipping sample %d: no sweep .mat files in %s\n', sample_id, folder_name);
    continue
end
fileT = nan(numel(mat_files),1); fileW = fileT; fileP = fileT;
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

% snap each requested value to the nearest available, then pick that file
uT = unique(fileT);  Tval = uT(find(abs(uT - T_target)  == min(abs(uT - T_target)),  1));
uW = unique(fileW);  Wval = uW(find(abs(uW - tw_target) == min(abs(uW - tw_target)), 1));
uP = unique(fileP);  Pval = uP(find(abs(uP - p_target)  == min(abs(uP - p_target)),  1));
if abs(Tval - T_target)  > 1e-9, fprintf('T:  requested %g -> %g\n',  T_target,  Tval); end
if abs(Wval - tw_target) > 1e-9, fprintf('tw: requested %g -> %g\n',  tw_target, Wval); end
if abs(Pval - p_target)  > 1e-9, fprintf('p:  requested %g -> %g\n',  p_target,  Pval); end

reltol = 1e-6;
match = find(abs(fileT - Tval) <= reltol*max(1,abs(Tval)) & ...
             abs(fileW - Wval) <= reltol*max(1,abs(Wval)) & ...
             abs(fileP - Pval) <= reltol*max(1,abs(Pval)), 1);
if isempty(match)
    fprintf('Skipping sample %d: no file for T=%g, tw=%g, p=%g\n', sample_id, Tval, Wval, Pval);
    continue
end
sim_data = load(fullfile(folder_name, mat_files(match).name));
xsaves = sim_data.xsaves;
Qsaves = sim_data.Qsaves;
num_tsim = size(xsaves, 1);
fprintf('==== sample %d (%s), t0=%d, %s/%s\n', sample_id, color_tag, t0, folder_name, mat_files(match).name);

% ---- which reference frames to show ---------------------------------------
if isempty(times_to_show)
    t_show = unique(round(linspace(t0, num_t, n_times)));
else
    t_show = unique(round(times_to_show));
end
t_show = t_show(t_show >= t0 & t_show <= num_t);
if isempty(t_show)
    fprintf('Skipping sample %d: no requested time in t0=%d..%d.\n', sample_id, t0, num_t);
    continue
end

% ---- interpolate the simulation onto those reference frames ---------------
% the n_saves snapshots span [0,T] <-> reference frames t0..num_t
sim_t = linspace(t0, num_t, num_tsim);
L_t0  = sum(sqrt(sum(diff(squeeze(data_x(t0,:,:)), 1, 1).^2, 2)));   % microns
xi = interp1(sim_t, squeeze(xsaves(:,:,1)), t_show) * L_t0;
yi = interp1(sim_t, squeeze(xsaves(:,:,2)), t_show) * L_t0;
zi = interp1(sim_t, squeeze(xsaves(:,:,3)), t_show) * L_t0;
d1xi = interp1(sim_t, squeeze(Qsaves(:,:,director_idx,1)), t_show);
d1yi = interp1(sim_t, squeeze(Qsaves(:,:,director_idx,2)), t_show);
d1zi = interp1(sim_t, squeeze(Qsaves(:,:,director_idx,3)), t_show);
% interp1 returns a row vector when a single time is requested
if numel(t_show) == 1
    xi = xi(:).';  yi = yi(:).';  zi = zi(:).';
    d1xi = d1xi(:).';  d1yi = d1yi(:).';  d1zi = d1zi(:).';
end

% ---- figure: rows = views, cols = time ------------------------------------
n_rows = numel(views_to_plot);
n_cols = numel(t_show);
fig = figure('Position', [10, 10, 240*n_cols+120, 210*n_rows+140]);
tiledlayout(n_rows, n_cols, 'Padding', 'tight', 'TileSpacing', 'compact');

for iv = 1:n_rows
    this_view = views_to_plot{iv};
    for it = 1:n_cols
        t = t_show(it);
        nexttile;

        trisurf(triangulation(meshes{t}.f, meshes{t}.v), 'edgecolor', 'none', ...
                'facecolor', 0.9*[1,1,1], 'faceAlpha', 0.2);
        hold on

        % imaged rod at this frame
        [xx, yy, zz] = tubeplot(squeeze(data_x(t,:,:))', tube_r, tube_n);
        surf(xx, yy, zz, 'Facecolor', 0.6*[1,1,1], 'EdgeColor', 'none');

        % simulation at this frame
        sim_x = [xi(it,:); yi(it,:); zi(it,:)];
        [xx, yy, zz] = tubeplot(sim_x, tube_r, tube_n);
        surf(xx, yy, zz, 'Facecolor', vik_c, 'EdgeColor', 'none');

        if draw_quivers
            % Every quiver3 argument is forced to a COLUMN with (:). MATLAB
            % drops trailing singleton dimensions, so data_x(t,rsub,1) comes
            % back as a 1xN row while d1r(rsub,1) is an Nx1 column, and mixing
            % the two gives "Z and U must be the same size".
            %
            % imaged frame: bases are element start nodes, so drop the last node
            d1r  = squeeze(data_Q(t, :, director_idx, :));   % (N-1) x 3
            rsub = 1:quiver_stride:size(d1r, 1);
            if draw_ref_frame
                rbx = data_x(t, rsub, 1);  rby = data_x(t, rsub, 2);
                rbz = data_x(t, rsub, 3);
                quiver3(rbx(:), rby(:), rbz(:), ...
                        d1r(rsub,1)*quiver_len, d1r(rsub,2)*quiver_len, ...
                        d1r(rsub,3)*quiver_len, ...
                        0, 'color', 'k', 'linewidth', quiver_lw);
            end
            % simulated frame
            ssub = 1:quiver_stride:size(d1xi, 2);
            sbx = xi(it,ssub);  sby = yi(it,ssub);  sbz = zi(it,ssub);
            sux = d1xi(it,ssub); suy = d1yi(it,ssub); suz = d1zi(it,ssub);
            quiver3(sbx(:), sby(:), sbz(:), ...
                    sux(:)*quiver_len, suy(:)*quiver_len, suz(:)*quiver_len, ...
                    0, 'color', vik_c, 'linewidth', quiver_lw);
            % the trailing 0 disables per-axes autoscaling, so quiver_len is
            % honoured identically in every tile
        end

        axis([xmin, xmax, ymin, ymax, zmin, zmax])
        daspect([1,1,1])
        axis off
        title(sprintf('t = %g min\n%s', (t - t0)*dt_min, this_view), 'fontsize', 7);

        lighting gouraud
        material dull
        if strcmp(this_view, 'side')
            view([0, 0]);
        elseif strcmp(this_view, 'dorsal')
            view(0, 90);
        else
            view(0, -90);
        end
        camlight('headlight');
    end
end

if draw_quivers
    qtag = sprintf('d%d frame: black = imaged, colour = simulation', director_idx);
else
    qtag = 'tubes only';
end
sgtitle(sprintf(['sample %d (%s), t0=%d,  T=%g, B_t/B=%g, Pl=%g' ...
                 '   (grey = imaged, colour = simulation)\n%s'], ...
                sample_id, color_tag, t0, Tval, Wval, Pval, qtag), 'fontsize', 9);

save_filename = sprintf('%s_dynamics_T%g_tw%g_p%g.png', char(folder_name), Tval, Wval, Pval);
print(save_filename, '-dpng', save_dpi);
fprintf('Saved %s\n', save_filename);
close(fig)

end   % samples_to_plot

fprintf('~done~\n');
