% fig1_coordinates_movie.m
% Movie frames of the surface coordinates on the gut mesh over time: left,
% the longitudinal coordinate (mesh hoop index, A -> P) in copper; right, the
% azimuthal coordinate (position around the hoop) in sky. Myocyte stripes in
% yellow, mesh end caps in grey, 50 um scale bar in the right panel. Fixed
% camera and limits, so the scale is the same in every frame. One PNG per
% timepoint in out_dir.

clear; close all
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..')); P = config_paths();

%%% settings
idx_sample = 4;           % WT #2
out_dir    = fullfile(P.out_fig1, sprintf('coords_s%d', idx_sample));
cam_view   = [-37.5, 30];
fact       = 1.052;       % push stripes slightly off the surface
stripe_w   = 6;           % myocyte stripe line width [pt]
cap_w      = 4;           % end-cap line width [pt]

%%% data
all_data = load(P.framesV6).all_data;
data = all_data{idx_sample};
num_t = min([numel(data.meshes), size(data.cent,1), size(data.myocytes_raw.red,1)]);
time_lags = get_time_lags(P, all_data);
lag = time_lags(idx_sample);

%%% fixed limits from every mesh
all_v = cell2mat(cellfun(@(m) m.v, data.meshes(1:num_t), 'UniformOutput', false));   % stacked vertices, all frames
lims = [min(all_v) - 5; max(all_v) + 5];

%%% panel layout and a fixed camera framing the gut tightly
pos   = {[0.005, 0.0, 0.49, 0.85], [0.505, 0.0, 0.49, 0.85]};   % normalized
fig_px = [1600, 650];
pad_frac = 0.04;          % margin around the gut, fraction of panel height

[V, U] = meshgrid(1:100, 1:100);   % U: hoop (longitudinal), V: around (azimuthal)
col_stripe = [0.929 0.694 0.125];
cmaps  = {copper, sky};
coords = {U, V};
titles = {'Longitudinal coordinate', 'Azimuthal coordinate'};

if ~exist(out_dir, 'dir'); mkdir(out_dir); end
fig = figure('Color', 'w', 'Position', [100, 100, fig_px]);

% camera: orthographic view(cam_view); target and view angle chosen so the
% projection of all meshes (all frames) just fills the panel
ax0 = axes(fig, 'Position', pos{1});
daspect(ax0, [1 1 1]); xlim(ax0, lims(:,1)); ylim(ax0, lims(:,2)); zlim(ax0, lims(:,3));
view(ax0, cam_view)
cpos = campos(ax0); ctar = camtarget(ax0); cup = camup(ax0);
fwd = (ctar - cpos) / norm(ctar - cpos);
rgt = cross(fwd, cup); rgt = rgt / norm(rgt);
upv = cross(rgt, fwd);
pr = (all_v - ctar) * rgt';  pu = (all_v - ctar) * upv';
c_r = 0.5*(min(pr) + max(pr)); c_u = 0.5*(min(pu) + max(pu));
half_w = 0.5*(max(pr) - min(pr)); half_h = 0.5*(max(pu) - min(pu));
panel_aspect = (pos{1}(3)*fig_px(1)) / (pos{1}(4)*fig_px(2));     % width / height
half_view = (1 + pad_frac) * max(half_h, half_w / panel_aspect);
cam_target = ctar + c_r*rgt + c_u*upv;
cam_pos    = cpos + c_r*rgt + c_u*upv;
cam_va     = 2*atand(half_view / norm(ctar - cpos));

%%% 50 um scale bar, in the screen plane (true length in the orthographic
%%% view), at the lower right of the gut's extent over all frames, moved
%%% towards the camera so the mesh never hides it
bar_len = 50;    % [um]
bar_r = max(pr) - [bar_len; 0];
bar_u = min(pu) + 0.04 * (max(pu) - min(pu));
bar_xyz = ctar + bar_r * rgt + bar_u * upv - 0.5 * norm(ctar - cpos) * fwd;   % 2 x 3
delete(ax0)

for t = 1:num_t
    clf
    msh = data.meshes{t};
    vert_UV = reshape(msh.v, [100,100,3]);
    cent = squeeze(data.cent(t,:,:));
    stripes = {data.myocytes_raw.red, data.myocytes_raw.green, ...
               data.myocytes_raw.yellow, data.myocytes_raw.blue};
    for p = 1:2
        ax = axes(fig, 'Position', pos{p});
        trisurf(triangulation(msh.f, msh.v), coords{p}(:), 'EdgeColor', 'none', 'FaceAlpha', 0.95);
        colormap(ax, cmaps{p})
        hold on
        for j = 1:4
            c = squeeze(stripes{j}(t,:,:));
            sp = cent + fact*(c - cent);
            plot3(sp(:,1), sp(:,2), sp(:,3), 'LineWidth', stripe_w, 'Color', col_stripe)
        end
        for hoop = [1, 100]   % end caps
            plot3(vert_UV(hoop,:,1), vert_UV(hoop,:,2), vert_UV(hoop,:,3), 'LineWidth', cap_w, 'Color', 0.5*[1,1,1])
        end
        if p == 2, plot3(bar_xyz(:,1), bar_xyz(:,2), bar_xyz(:,3), 'k', 'LineWidth', 3, 'Clipping', 'off'), end
        lighting gouraud; material dull
        daspect([1,1,1])
        xlim(lims(:,1)); ylim(lims(:,2)); zlim(lims(:,3))
        view(cam_view)
        set(ax, 'CameraTarget', cam_target, 'CameraPosition', cam_pos, ...
                'CameraUpVector', cup, 'CameraViewAngle', cam_va)
        axis off
        % label at a fixed figure position (axes titles follow the manual camera)
        annotation(fig, 'textbox', [pos{p}(1), 0.86, pos{p}(3), 0.06], 'String', titles{p}, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14)
    end
    sgtitle(sprintf('%d min', round(aligned_time(t, lag, 'frame', P))))
    save_frame(fig, fullfile(out_dir, sprintf('%06d.png', t)), 200)
end
