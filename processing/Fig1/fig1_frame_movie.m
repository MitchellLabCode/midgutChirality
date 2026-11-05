% fig1_frame_movie.m
% Movie frames of the material frame along the centerline over time: the
% gut mesh (transparent), the centerline as a tube, and the frame vectors
% d1 (L-R) and d2 (V-D) at every centerline edge. Fixed camera, limits and
% arrow length, so the scale is the same in every frame. One PNG per
% timepoint in out_dir.

clear; close all
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..')); P = config_paths();

%%% settings
idx_sample = 4;           % WT #2
out_dir    = fullfile(P.out_fig1, sprintf('frame_s%d', idx_sample));
cam_view   = [-37.5, 30];
arrow_len  = 18;          % frame-vector length [um] (fixed, no quiver autoscale)
tube_r     = 1.5;         % centerline tube radius [um]
show_d3    = false;       % also draw the tangent d3
cl = lines(4);
dcol = {cl(2,:), cl(3,:), cl(1,:)};    % d1 red-orange, d2 yellow, d3 blue

%%% data
all_data = load(P.framesV6).all_data;
data = all_data{idx_sample};
num_t = min([numel(data.meshes), size(data.cent,1), size(data.Q,1)]);
time_lags = get_time_lags(P, all_data);
lag = time_lags(idx_sample);

all_v = cell2mat(cellfun(@(m) m.v, data.meshes(1:num_t), 'UniformOutput', false));   % stacked vertices, all frames
lims = [min(all_v) - 5; max(all_v) + 5];
dirs = [1 2]; if show_d3, dirs = [1 2 3]; end

if ~exist(out_dir, 'dir'); mkdir(out_dir); end
fig = figure('Color', 'w', 'Position', [100, 100, 1300, 820]);
for t = 1:num_t
    clf
    msh  = data.meshes{t};
    cent = squeeze(data.cent(t,:,:));
    Q    = squeeze(data.Q(t,:,:,:));             % (N-1) x 3 x 3, Q(:,i,:) = d_i
    base = cent(1:end-1,:);                      % frame lives on edges; drawn at the edge start

    trisurf(triangulation(msh.f, msh.v), 'EdgeColor', 'none', 'FaceColor', 0.9*[1,1,1], 'FaceAlpha', 0.1);
    hold on
    [tx, ty, tz] = tubeplot(cent', tube_r, 16, 1e-9);
    surf(tx, ty, tz, 'FaceColor', 0.3*[1,1,1], 'EdgeColor', 'none')
    for i = dirs
        d = squeeze(Q(:,i,:));
        quiver3(base(:,1), base(:,2), base(:,3), arrow_len*d(:,1), arrow_len*d(:,2), arrow_len*d(:,3), 0, ...
            'LineWidth', 2, 'Color', dcol{i}, 'MaxHeadSize', 0.5)
    end
    daspect([1,1,1])
    xlim(lims(:,1)); ylim(lims(:,2)); zlim(lims(:,3))
    view(cam_view)
    axis off
    lighting gouraud; material dull; camlight
    title(sprintf('%d min', round(aligned_time(t, lag, 'frame', P))), 'FontSize', 16, 'FontWeight', 'normal')
    save_frame(fig, fullfile(out_dir, sprintf('%06d.png', t)), 200)
end
