% fig1_chiral_movie.m
% Movie frames of the chirality q on the gut: a WT and a byn>Myo1C embryo
% side by side (columns), ventral and side views (rows), on the aligned time
% axis. Same measure and colours as fig1_chirality_10samples.m (q over 2w+1
% tangents, Savitzky-Golay smoothed along s, managua_br, limits +-1), the
% centerline drawn as a tube coloured by q. Frames run over the aligned
% times both embryos cover; one PNG per frame in out_dir. Assemble with
% make_movie.py.

clear; close all
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..')); P = config_paths();

%%% which embryos and frames
% (the Fig. 1 movie, see make_movie.py: WT_idx = 4 and
%  fig_title = 'The resulting shapes display \itperverted\rm chirality')
WT_idx  = 6;     % WT #3
Myo_idx = 1;     % byn>Myo1C #1
panel_titles = {sample_label(P, WT_idx), sample_label(P, Myo_idx)};   % column titles
out_dir = fullfile(P.out_fig1, sprintf('chiral3d_s%d_s%d', WT_idx, Myo_idx));
which_frames = 'all';   % 'all', 'last' (last time both embryos share), or aligned times [min]
fig_title = '';         % optional title above the timestamp (TeX, e.g. '\itword\rm'); '' for none

all_data = load(P.framesV6).all_data;
time_lags = round(get_time_lags(P, all_data));   % frame t of sample s is aligned frame t + lag(s)

samples_indices = [WT_idx, Myo_idx];
w = 5;            % q window half-width
tube_r = 2.5;     % centerline tube radius [um]
cmap = managua_br();
cl = [-1, 1];

x_interp = cell(1, 2); num_t = zeros(1, 2);
for p = 1:2
    x_interp{p} = arclength_resample(all_data{samples_indices(p)}.x);
    num_t(p) = size(x_interp{p}, 1);
end

% aligned frames covered by both embryos
lags = time_lags(samples_indices);
a_frames = max(1 + lags) : min(num_t + lags);
a_first  = a_frames(1);           % file names index the full sequence
if ischar(which_frames) && strcmp(which_frames, 'last')
    a_frames = a_frames(end);
elseif isnumeric(which_frames)
    a_frames = a_frames(ismember(aligned_time(a_frames, 0, 'frame', P), which_frames));
end
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

% views; anterior (low x) on the left in both
view_names = {'ventral', 'side'};
view_angles = {[0, -90], [0, 0]};
% panel layout (normalized figure units): rows = views, columns = embryos
col_x = [0.07, 0.49];  col_w = 0.40;
row_y = [0.50, 0.08];  row_h = 0.38;
if ~isempty(fig_title), row_y = [0.46, 0.06]; row_h = 0.36; end

fig = figure('Color', 'w', 'Position', [100, 100, 1400, 900]);
for k = 1:numel(a_frames)
    a = a_frames(k);
    t_min = aligned_time(a, 0, 'frame', P);           % = aligned_time(tidx, lag, 'frame')
    clf
    for p = 1:2
        data = all_data{samples_indices(p)};
        tidx = a - lags(p);
        msh  = data.meshes{tidx};
        curv = squeeze(x_interp{p}(tidx, 2:end, :));
        [q, centres] = coiling_profile(curv, w);
        q = savgol(q, 2, 5);

        for r = 1:2
            ax = axes(fig, 'Position', [col_x(p), row_y(r), col_w, row_h]);
            trisurf(triangulation(msh.f, msh.v), 'EdgeColor', 'none', ...
                'FaceColor', 0.9*[1,1,1], 'FaceAlpha', 0.35);
            hold on
            % tubeplot adds one ring at each end, hence the padded colours
            [tx, ty, tz] = tubeplot(curv(centres,:)', tube_r, 16, 1e-9);
            tc = repmat([q(1); q(:); q(end)]', size(tx,1), 1);
            surf(tx, ty, tz, tc, 'EdgeColor', 'none', 'FaceColor', 'interp')
            colormap(ax, cmap); clim(ax, cl)
            axis equal off
            xlim([-50, 275]); ylim([-85, 85]); zlim([-75, 80])
            lighting gouraud; material dull
            view(view_angles{r})
            camlight
            if r == 1
                camlight(0, -90)
                title(panel_titles{p}, 'FontSize', 14)
            end
            if r == 2 && p == 2      % 50 um scale bar along x
                plot3([200, 250], -85*[1, 1], -70*[1, 1], 'Color', 'k', 'LineWidth', 2)
            end
        end
    end
    for r = 1:2
        annotation(fig, 'textbox', [0.0, row_y(r), 0.06, row_h], 'String', view_names{r}, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 13);
    end
    % colorbar with end caps: values beyond +-1 are shown saturated
    colorbar_extend(fig, [0.925, 0.33, 0.012, 0.34], cmap, cl, [-1, 0, 1], 'q');
    if isempty(fig_title)
        sgtitle(sprintf('%d min', t_min))
    else
        annotation(fig, 'textbox', [0, 0.92, 1, 0.07], 'String', fig_title, 'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 20)
        annotation(fig, 'textbox', [0, 0.865, 1, 0.05], 'String', sprintf('%d min', t_min), ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 15)
    end
    save_frame(fig, fullfile(out_dir, sprintf('%06d.png', a - a_first + 1)), 300)
end
