% fig1_tantrix_movie.m
% Movie frames illustrating the chirality measure. Left: the midgut (mesh
% and centerline tube) with the unit tangent moving along the centerline.
% Right: the same tangent drawn from the centre of a unit sphere, sweeping
% out the tangent indicatrix (tantrix). Both panels share the camera, so the
% arrows point the same way. One PNG per centerline position in out_dir.

clear; close all
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..')); P = config_paths();

%%% settings
idx_sample = 4;           % WT #2
t_frame    = 60;          % time frame
out_dir    = fullfile(P.out_fig1, sprintf('tantrix_s%d', idx_sample));
cam_view   = [-37.5, 30];
arrow_len  = 40;                 % tangent arrow length on the gut [um]
col_cent   = 0.25*[1 1 1];       % centerline tube (dark grey)
col_t      = [0 0.447 0.741];    % tangent vector (blue)
col_trace  = [0 0 0];            % tantrix on the sphere (black)
fig_title  = 'We define a chirality measure through the solid angle swept by the tangent vector';
trace_r = 0.012;      % tantrix tube radius (unit sphere)
tube_r  = 1.5;        % centerline tube radius [um]
base_r  = 2.2;        % sphere marking the tangent's base [um]

%%% data
all_data = load(P.framesV6).all_data;
data = all_data{idx_sample};
msh  = data.meshes{t_frame};
cent = squeeze(data.cent(t_frame, :, :));
tang = diff(cent, 1, 1);
tang = tang ./ vecnorm(tang, 2, 2);                 % unit tangents
mid  = 0.5 * (cent(1:end-1,:) + cent(2:end,:));     % tangent base points
n_tan = size(tang, 1);

% fixed limits for the gut panel, covering the mesh and every arrow tip
pts = [msh.v; mid + arrow_len*tang];
pad = 5;   % um
lim_gut = [min(pts) - pad; max(pts) + pad];

if ~exist(out_dir, 'dir'); mkdir(out_dir); end
[sx, sy, sz] = sphere(64);
th = linspace(0, 2*pi, 200);

fig = figure('Color', 'w', 'Position', [100, 100, 1600, 800]);
for k = 1:n_tan
    clf
    %%% left: midgut with the moving tangent
    axes(fig, 'Position', [0.02, 0.05, 0.52, 0.80]);
    trisurf(triangulation(msh.f, msh.v), 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'FaceColor', 0.9*[1,1,1]);
    hold on
    [tx, ty, tz] = tubeplot(cent', tube_r, 16, 1e-9);
    surf(tx, ty, tz, 'FaceColor', col_cent, 'EdgeColor', 'none')
    arrow3d(mid(k,:), arrow_len*tang(k,:), 1.2, 3.2, 10, col_t)
    surf(base_r*sx + mid(k,1), base_r*sy + mid(k,2), base_r*sz + mid(k,3), 'FaceColor', col_t, 'EdgeColor', 'none')
    axis equal off
    xlim(lim_gut(:,1)); ylim(lim_gut(:,2)); zlim(lim_gut(:,3))
    view(cam_view)
    lighting gouraud; material dull; camlight

    %%% right: unit sphere with the tangent from the origin and its trace
    axes(fig, 'Position', [0.56, 0.08, 0.40, 0.74]);
    surf(sx, sy, sz, 'FaceColor', 0.85*[1,1,1], 'FaceAlpha', 0.3, 'EdgeColor', 'none')
    hold on
    gc = 0.75*[1 1 1];                                % equator and meridians as depth cues
    plot3(cos(th), sin(th), 0*th, 'Color', gc, 'LineWidth', 0.5)
    plot3(cos(th), 0*th, sin(th), 'Color', gc, 'LineWidth', 0.5)
    plot3(0*th, cos(th), sin(th), 'Color', gc, 'LineWidth', 0.5)
    % tantrix so far, as a thin lit tube just outside the sphere, so it is
    % shaded like the sphere and fades behind it on the far hemisphere
    if k >= 2
        [qx, qy, qz] = tubeplot(1.015 * tang(1:k,:)', trace_r, 12, 1e-9);
        surf(qx, qy, qz, 'FaceColor', col_trace, 'EdgeColor', 'none')
    end
    arrow3d([0 0 0], tang(k,:), 0.025, 0.07, 0.2, col_t)
    axis equal off
    xlim([-1.05 1.05]); ylim([-1.05 1.05]); zlim([-1.05 1.05])
    view(cam_view)
    lighting gouraud; material dull; camlight

    annotation(fig, 'textbox', [0, 0.88, 1, 0.1], 'String', fig_title, 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 18)
    save_frame(fig, fullfile(out_dir, sprintf('%06d.png', k)), 200)
end

