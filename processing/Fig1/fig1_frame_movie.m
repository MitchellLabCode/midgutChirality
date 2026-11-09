% fig1_frame_movie.m
% Movie frames of the material frame along the centerline over time: the
% gut mesh (transparent), the centerline as a tube, and the frame vectors
% d1 (L-R) and d2 (V-D) at every centerline edge. Fixed camera, limits and
% arrow length, so the scale is the same in every frame; 50 um scale bar.
% One PNG per timepoint in out_dir; timestamp top left, minutes from the
% first frame.

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

all_v = cell2mat(cellfun(@(m) m.v, data.meshes(1:num_t), 'UniformOutput', false));   % stacked vertices, all frames
lims = [min(all_v) - 5; max(all_v) + 5];
dirs = [1 2]; if show_d3, dirs = [1 2 3]; end

if ~exist(out_dir, 'dir'); mkdir(out_dir); end
ts_font = 'Helvetica Neue';  % timestamp font; falls back to Nimbus Sans (Helvetica clone) if not installed
if ~any(strcmp(listfonts, ts_font)), ts_font = 'Nimbus Sans'; end
fig = figure('Color', 'w', 'Position', [100, 100, 1300, 820]);

% 50 um scale bar, in the screen plane (true length in the orthographic
% view), at the lower right of the gut's extent over all frames, moved
% towards the camera so the mesh never hides it. The camera is the one the
% loop below gets from view(cam_view) with these limits.
ax0 = axes(fig);
daspect(ax0, [1 1 1]); xlim(ax0, lims(:,1)); ylim(ax0, lims(:,2)); zlim(ax0, lims(:,3));
view(ax0, cam_view)
cpos = campos(ax0); ctar = camtarget(ax0); cup = camup(ax0);
delete(ax0)
fwd = (ctar - cpos) / norm(ctar - cpos);
rgt = cross(fwd, cup); rgt = rgt / norm(rgt);
upv = cross(rgt, fwd);
pr = (all_v - ctar) * rgt';  pu = (all_v - ctar) * upv';
bar_len = 50;    % [um]
bar_r = max(pr) - [bar_len; 0];
bar_u = min(pu) + 0.04 * (max(pu) - min(pu));
bar_xyz = ctar + bar_r * rgt + bar_u * upv - 0.5 * norm(ctar - cpos) * fwd;   % 2 x 3
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
    plot3(bar_xyz(:,1), bar_xyz(:,2), bar_xyz(:,3), 'k', 'LineWidth', 3, 'Clipping', 'off')
    daspect([1,1,1])
    xlim(lims(:,1)); ylim(lims(:,2)); zlim(lims(:,3))
    view(cam_view)
    axis off
    lighting gouraud; material dull; camlight
    annotation(fig, 'textbox', [0.03, 0.84, 0.2, 0.06], 'String', sprintf('%d min', P.dt*(t-1)), ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'FontName', ts_font, 'FontSize', 16)   % timestamp, top left, from t = 0
    save_frame(fig, fullfile(out_dir, sprintf('%06d.png', t)), 200)
end
