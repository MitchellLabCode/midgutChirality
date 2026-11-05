% asym_constrictions.m
% Gut surface of each embryo at given aligned times, with the three
% constrictions (fold hoops, coloured, with their centroids) and the hoop
% centroid line (grey); ventral and side views, all panels rendered in one
% shared box so that sizes are comparable across embryos and views.
% Writes <WT|Myo>_asym_<ventral|side>_small_s<N>_t<T>.png to P.out_fig2/asym.
%
% Frame selection: with P.time_convention = 'published' this reproduces the
% published panels, whose frame lookup also added the 2nd-fold onset of the
% reference embryo (fold_onset(2) of sample 4 = 9 frames) on top of the lags;
% the panel labelled T min is then at T + 18 min on the axis of the other
% figures. With 'unified' the frame at aligned time T (aligned_time 'frame')
% is used.

clear; close all
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..')); P = config_paths();
outdir = fullfile(P.out_fig2, 'asym');
if ~exist(outdir, 'dir'); mkdir(outdir); end

all_data = load(P.framesV6).all_data;
time_lags = get_time_lags(P, all_data);

%%% parameters
samples_to_plot = 1:10;
aligned_times   = [20, 30];       % [min]
if strcmp(P.time_convention, 'published')
    frame_offset = all_data{4}.fold_onset(2);
else
    frame_offset = 0;
end

%%% shared axis limits: every panel is rendered in the same box
xl = [-50, 275]; yl = [-85, 85]; zl = [-75, 90];

%%% export settings (as the interactive "Export To..." dialog)
figSize_px  = [560, 420];
exportWidth = 2.1;          % inches
exportDPI   = 600;

frame_at = @(T, i) max(1, min(numel(all_data{i}.meshes), ...
    round(T / P.dt - time_lags(i)) + frame_offset));   % inverse of aligned_time 'frame'

%%% check the shared box against every mesh that will be drawn
bbox = [inf, -inf; inf, -inf; inf, -inf];
for i = samples_to_plot
    for T = aligned_times
        v = all_data{i}.meshes{frame_at(T, i)}.v;
        bbox(:,1) = min(bbox(:,1), min(v, [], 1)');
        bbox(:,2) = max(bbox(:,2), max(v, [], 1)');
    end
end
clipped = [xl; yl; zl] - bbox;
if any(clipped(:,1) > 0) || any(clipped(:,2) < 0)
    warning('asym:clippedLimits', 'Shared limits clip the data: x [%.4g %.4g], y [%.4g %.4g], z [%.4g %.4g]', ...
        bbox(1,:), bbox(2,:), bbox(3,:));
end

cmap_lines = lines(3);
for i = samples_to_plot
    data = all_data{i};
    if P.samples(i).is_WT, geno = 'WT'; else, geno = 'Myo'; end
    for T = aligned_times
        t = frame_at(T, i);
        msh = data.meshes{t};
        vert_UV = reshape(msh.v, [100, 100, 3]);

        figure('Units', 'pixels', 'Position', [100, 100, figSize_px], 'Color', 'w')
        trisurf(triangulation(msh.f, msh.v), 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'FaceColor', 0.9*[1,1,1])
        hold on
        title(sprintf('Sample %d | aligned time = %.0f min', i, T));
        for k = 1:3                                     % constrictions
            hoop = squeeze(vert_UV(data.folds(t, k), :, :));
            plot3(hoop(:,1), hoop(:,2), hoop(:,3), 'LineWidth', 5, 'Color', cmap_lines(k,:))
            scatter3(mean(hoop(:,1)), mean(hoop(:,2)), mean(hoop(:,3)), 'filled', 'MarkerFaceColor', cmap_lines(k,:))
        end
        plot3(mean(vert_UV(:,:,1),2), mean(vert_UV(:,:,2),2), mean(vert_UV(:,:,3),2), ...
            'Color', 0.7*[1,1,1], 'LineWidth', 4)
        lighting gouraud; material dull; axis off

        h = camlight(0, -90);
        view(0, -90)
        applySharedAxes(gca, xl, yl, zl)
        exportgraphics(gcf, fullfile(outdir, sprintf('%s_asym_ventral_small_s%d_t%d.png', geno, i, T)), ...
            'Width', exportWidth, 'Height', 'auto', 'Units', 'inches', 'Resolution', exportDPI, 'Padding', 'figure')
        delete(h); camlight
        view(0, 0)
        applySharedAxes(gca, xl, yl, zl)
        exportgraphics(gcf, fullfile(outdir, sprintf('%s_asym_side_small_s%d_t%d.png', geno, i, T)), ...
            'Width', exportWidth, 'Height', 'auto', 'Units', 'inches', 'Resolution', exportDPI, 'Padding', 'figure')
        close(gcf)
    end
end


function applySharedAxes(ax, xl, yl, zl)
% pin limits, aspect and zoom so every panel renders at the same scale
set(ax, 'XLim', xl, 'YLim', yl, 'ZLim', zl, 'DataAspectRatio', [1, 1, 1], ...
    'PlotBoxAspectRatio', [diff(xl), diff(yl), diff(zl)], 'CameraViewAngleMode', 'auto');
set(ax, 'CameraViewAngleMode', 'manual');
end
