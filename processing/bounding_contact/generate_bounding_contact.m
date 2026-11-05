% generate_bounding_contact.m
% Reference data and confinement for the Cosserat-rod simulations (cosserat/):
% for every embryo, the bounding walls and the contact radius of the gut.
%
%   walls          : the convex hull of every second centerline point over all
%                    times, as half-spaces n . x < d (boundingplanes.m)
%   contact radius : 1.2 x the distance from each centerline point to the
%                    nearest surface vertex, at every time (contact_radius.m)
%
% Writes P.bounding_contact/refIV6_sample<N>_walls3_cball.mat (cosserat/data/)
% with x (centerline), Q (frame), radii (hoop radii), n, d (walls), tlag
% (time lag [frames]), fold3 (mesh-hoop index of the 3rd fold per frame, -1
% before its onset), circum_radii (contact radius), and the contact-radius
% heatmaps to P.out_fig2.

clear; close all
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'), here); P = config_paths();

save_refs  = true;     % write the reference files
show_plots = false;    % per-sample 3D views of the walls and contact spheres

all_data = load(P.framesV6).all_data;
time_lags = get_time_lags(P, all_data);
if ~exist(P.bounding_contact, 'dir'); mkdir(P.bounding_contact); end
if ~exist(P.out_fig2, 'dir'); mkdir(P.out_fig2); end

for idx_sample = 1:numel(all_data)
    cdat = all_data{idx_sample};
    x = cdat.cent;
    Q = cdat.Q;
    radii = cdat.radii;
    fold3 = double(cdat.folds(:,3)); fold3(1:cdat.fold_onset(3)) = -1;
    tlag = time_lags(idx_sample);

    %%% walls: convex hull of every 2nd centerline point, all times
    dx = x(:,1:2:end,1); dy = x(:,1:2:end,2); dz = x(:,1:2:end,3);
    [n, d] = boundingplanes(dx(:), dy(:), dz(:));

    %%% contact radius
    circum_radii = contact_radius(cdat);
    [num_t, num_s, ~] = size(x);

    if show_plots
        figure('Color', 'w')
        k2 = convhull(dx(:), dy(:), dz(:), 'Simplify', true);
        scatter3(dx(:), dy(:), dz(:)); hold on
        trisurf(k2, dx(:), dy(:), dz(:), 'FaceColor', 'cyan')
        verts = [dx(:), dy(:), dz(:)];
        fc = (verts(k2(:,1),:) + verts(k2(:,2),:) + verts(k2(:,3),:)) / 3;
        quiver3(fc(:,1), fc(:,2), fc(:,3), n(:,1), n(:,2), n(:,3))
        daspect([1,1,1])

        figure('Color', 'w')
        t0 = num_t;
        plot3(x(t0,:,1), x(t0,:,2), x(t0,:,3)); hold on
        trisurf(triangulation(cdat.meshes{t0}.f, cdat.meshes{t0}.v), 'edgecolor', 'none', ...
            'facecolor', 0.9*[1,1,1], 'faceAlpha', 0.3);
        [xs, ys, zs] = sphere(30);
        for s = 2:num_s
            surf(circum_radii(t0,s)*xs + x(t0,s,1), circum_radii(t0,s)*ys + x(t0,s,2), ...
                 circum_radii(t0,s)*zs + x(t0,s,3), 'FaceColor', 'r', 'edgecolor', 'none')
        end
        daspect([1,1,1]); view(0,-90); camlight('headlight')
        lighting gouraud; material dull; axis off
    end

    %%% heatmap of the contact radius, with folds and the A/P split
    figure('Color', 'w')
    colormap(bone)
    t_al = aligned_time((1:num_t)', tlag, 'frame', P);
    imagesc('XData', linspace(0, 1, num_s), 'YData', t_al, 'CData', circum_radii)
    hold on
    folds = cdat.folds;
    for k = 1:3, folds(1:cdat.fold_onset(k), k) = -1; end
    fold_s = fold_arclength(cdat, 'cent');
    for k = 1:3
        valid = folds(:,k) > 0;
        plot(fold_s(valid,k), t_al(valid), 'Color', 0.6*[1,1,1], 'LineWidth', 1)
    end
    idx_split = zeros(num_t, 1);
    for t = 1:num_t, [~, idx_split(t)] = max(x(t,:,1)); end
    plot((idx_split - 1)/(num_s - 1), t_al, 'Color', 0.6*[1,1,1], 'LineStyle', '--')
    xlim([0, 1]); ylim([t_al(1), t_al(end)])
    colorbar; clim([0, 75])
    set(gca, 'YDir', 'reverse')
    title(['Sample ' num2str(idx_sample)])
    exportgraphics(gcf, fullfile(P.out_fig2, sprintf('cradius_heatmap_sample%d.pdf', idx_sample)), 'ContentType', 'vector')

    if save_refs
        save(fullfile(P.bounding_contact, sprintf('refIV6_sample%d_walls3_cball.mat', idx_sample)), ...
            'x', 'Q', 'radii', 'n', 'd', 'tlag', 'fold3', 'circum_radii')
    end
    fprintf('sample %d done\n', idx_sample);
end
