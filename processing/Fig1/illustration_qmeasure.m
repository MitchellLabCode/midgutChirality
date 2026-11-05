% illustration_qmeasure.m
% Illustration of the chirality measure q (Fig. 1): left, a stretch of gut
% centerline with three unit tangents; right, the tangents on the unit
% sphere (tantrix) and the signed spherical area q_h enclosed by 2w+1
% consecutive tangents (qmeasure.m), shaded in red.
% Writes illustration_qmeasure_curve.pdf and illustration_qmeasure_sphere.pdf
% to P.out_fig1. The annotated movie version is fig1_qmeasure_frame.m.

clear; close all
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..')); P = config_paths();
if ~exist(P.out_fig1, 'dir'); mkdir(P.out_fig1); end

tidx = 50;                       % time frame of the example embryo
i_curve  = [55, 65, 80];         % tangents drawn on the curve
i_sphere = [60, 65, 80];         % tangents drawn on the sphere
i_patch = 70; half_w = 7;        % centre and half-width of the shaded q_h window
cols = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]};

%%% centerline and unit tangents
dat = load(P.qmeasure_example, 'dat').dat;
curv = dat{tidx}.mesh.avgpts_um(10:end-1,:);
curv = smooth_curve3D(curv, 1, size(curv,1), 1e-3, 3);
tangents = diff(curv, [], 1);
tangents = tangents ./ vecnorm(tangents, 2, 2);

bp = [tangents(i_patch-half_w:i_patch+half_w,:); tangents(i_patch-half_w,:)];
fprintf('q_h at tangent %d (window %d): %.4f\n', i_patch, 2*half_w+1, ...
    qmeasure(tangents(i_patch-half_w:i_patch+half_w,:)));

%%% left: centerline as a tube with three tangents
figure('Color', 'w')
[x, y, z] = tubeplot(curv(50:end,:)', 1.5, 12);
surf(x, y, z, 'FaceColor', 0.9*[1,1,1], 'EdgeColor', 'none')
hold on
for k = 1:3
    i = i_curve(k);
    p0 = curv(i,:);
    if k == 2    % offset the middle arrow off the tube along the curvature direction
        nv = tangents(i+1,:) - tangents(i,:); nv = nv / norm(nv);
        p0 = p0 + 2.5 * [nv(1), -nv(2), nv(3)];
    end
    quiver3(p0(1), p0(2), p0(3), tangents(i,1), tangents(i,2), tangents(i,3), 20, ...
        'filled', 'LineWidth', 3, 'Color', cols{k})
end
daspect([1,1,1]); camlight; material dull; axis off
exportgraphics(gcf, fullfile(P.out_fig1, 'illustration_qmeasure_curve.pdf'), 'ContentType', 'vector')

%%% right: tantrix on the unit sphere and the q_h patch
figure('Color', 'w')
[xs, ys, zs] = sphere(500);
surf(xs, ys, zs, 'EdgeColor', 'none', 'FaceColor', 0.9*[1,1,1], 'FaceAlpha', 0.3)
hold on
plot3(tangents(50:end,1), tangents(50:end,2), tangents(50:end,3), 'LineWidth', 3, 'Color', 'k')
for k = 1:3
    i = i_sphere(k);
    quiver3(0, 0, 0, tangents(i,1), tangents(i,2), tangents(i,3), 1, 'filled', 'LineWidth', 3, 'Color', cols{k})
end
th = 0:0.01:2*pi;
plot3(cos(th), sin(th), zeros(size(th)), 'Color', 'k', 'LineWidth', 1)    % equator
[tri, Xq, Yq, Zq] = spherical_patch(bp);
trisurf(tri, 1.01*Xq, 1.01*Yq, 1.01*Zq, 'FaceColor', [1, 0, 0], 'FaceAlpha', 0.8, 'EdgeColor', 'none')
camlight; material dull
daspect([1,1,1])
view([135 90])
set(gca, 'Visible', 'off')
exportgraphics(gcf, fullfile(P.out_fig1, 'illustration_qmeasure_sphere.pdf'), 'ContentType', 'image', 'Resolution', 400)
