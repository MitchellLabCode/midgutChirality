% fig1_qmeasure_frame.m
% Still frame for the chirality movie (between the tantrix and the chiral
% movie parts, see make_movie.py): the tantrix of the example centerline on
% the unit sphere, three tangent vectors, and the signed spherical area
% q_h(s2), with labels and the movie title. Same data and construction as
% the sphere panel of illustration_qmeasure.m. Writes one PNG to out_dir.

clear; close all
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..')); P = config_paths();

%%% settings
tidx     = 50;
out_dir  = fullfile(P.out_fig1, 'qmeasure_frame');
cam_view = [141, 32];
title_str = 'The signed spherical area defines a measure of chirality';
col1 = [0, 0.4470, 0.7410];       % t(s1)
col2 = [0.8500, 0.3250, 0.0980];  % t(s2)
col3 = [0.9290, 0.6940, 0.1250];  % t(s3)
col_q = [0.65, 0.08, 0.10];       % q_h patch
i1 = 60; i2 = 65; i3 = 80;        % tangent indices, as in illustration_qmeasure.m
i_patch = 70; half_w = 7;         % patch window, as in illustration_qmeasure.m

%%% data and tangents
dat = load(P.qmeasure_example, 'dat').dat;
curv = dat{tidx}.mesh.avgpts_um(10:end-1,:);
curv = smooth_curve3D(curv, 1, size(curv,1), 1e-3, 3);
tangents = diff(curv, [], 1);
tangents = tangents ./ vecnorm(tangents, 2, 2);
tr = tangents(50:end, :);         % plotted part of the tantrix
bp = [tangents(i_patch-half_w:i_patch+half_w,:); tangents(i_patch-half_w,:)];
[tri, Xq, Yq, Zq] = spherical_patch(bp);

%%% figure
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
fig = figure('Color', 'w', 'Position', [100, 100, 1100, 1000]);
ax = axes(fig, 'Position', [0.05, 0.02, 0.9, 0.86]);
hold(ax, 'on')
[xs, ys, zs] = sphere(200);
surf(ax, xs, ys, zs, 'EdgeColor', 'none', 'FaceColor', 0.9*[1 1 1], 'FaceAlpha', 0.3)
th = linspace(0, 2*pi, 400);
plot3(ax, cos(th), sin(th), 0*th, 'k', 'LineWidth', 1)                       % equator
trisurf(tri, 1.01*Xq, 1.01*Yq, 1.01*Zq, 'FaceColor', col_q, 'FaceAlpha', 0.85, 'EdgeColor', 'none')
plot3(ax, 1.005*tr(:,1), 1.005*tr(:,2), 1.005*tr(:,3), 'k', 'LineWidth', 2)   % tantrix

% direction marks along the tantrix (small cones pointing along the curve)
for k = round(linspace(6, size(tr,1)-6, 4))
    d = tr(k+1,:) - tr(k-1,:);
    arrow3d(1.008*tr(k,:) - 0.05*d/norm(d), 0.1*d/norm(d), 0.0001, 0.035, 0.1, [0 0 0])
end
% the three tangent vectors from the origin
arrow3d([0 0 0], tangents(i1,:), 0.012, 0.04, 0.12, col1)
arrow3d([0 0 0], tangents(i2,:), 0.012, 0.04, 0.12, col2)
arrow3d([0 0 0], tangents(i3,:), 0.012, 0.04, 0.12, col3)
[sx, sy, sz] = sphere(30);
surf(ax, 0.025*sx, 0.025*sy, 0.025*sz, 'FaceColor', 'k', 'EdgeColor', 'none')  % origin

axis(ax, 'equal', 'off')
xlim(ax, [-1.1 1.1]); ylim(ax, [-1.1 1.1]); zlim(ax, [-1.1 1.1])
view(ax, cam_view)
camzoom(ax, 1.45)
lighting(ax, 'gouraud'); material(ax, 'dull'); camlight(ax)

% labels beside the vector tips / patch, offset in screen space (right, up)
v_r = [cosd(cam_view(1)), sind(cam_view(1)), 0];                     % screen right
v_c = [cosd(cam_view(2))*sind(cam_view(1)), -cosd(cam_view(2))*cosd(cam_view(1)), sind(cam_view(2))];
v_u = cross(v_r, v_c); v_u = v_u / norm(v_u);                         % screen up
if v_u(3) < 0, v_u = -v_u; end
lbl  = {'$\mathbf{t}(s_1)$', '$\mathbf{t}(s_2)$', '$q_h(s_2)$'};
pts  = {tangents(i1,:), tangents(i2,:), mean(bp(1:end-1,:))/norm(mean(bp(1:end-1,:)))};
offs = {[-0.12, -0.10], [-0.05, -0.17], [0.16, -0.02]};               % [right, up]
algn = {'right', 'center', 'left'};
cols = {col1, col2, col_q};
for k = 1:3
    q = pts{k} + offs{k}(1)*v_r + offs{k}(2)*v_u;
    text(ax, q(1), q(2), q(3), lbl{k}, 'Interpreter', 'latex', 'FontSize', 40, ...
        'Color', cols{k}, 'HorizontalAlignment', algn{k}, 'VerticalAlignment', 'middle')
end
annotation(fig, 'textbox', [0, 0.90, 1, 0.08], 'String', title_str, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 22)
save_frame(fig, fullfile(out_dir, '000001.png'), 200)
