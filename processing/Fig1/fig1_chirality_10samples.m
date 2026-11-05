% fig1_chirality_10samples.m
% Chirality of the midgut centerline for all 10 embryos (Fig. 1, Fig. S-chiral):
%
%   1. 3D renders at the first and last frame with a mesh, the centerline
%      coloured by the quasi-local chirality q (qmeasure.m over 2w+1 tangents):
%        chiral_sample<N>.pdf       q as dots
%        chiral_sample<N>_tube.png  q on a tube around the centerline
%   2. chirality kymographs q(s/L, t) on the aligned time axis, with the folds
%      and the anterior/posterior split (max-x point):    kymos_<N>.pdf
%   3. PCA of the pooled kymograph rows: pca_scatter.pdf, pca_loadings.pdf
%
% Outputs go to P.out_fig1/chirality. The kymographs are exported to the
% figshare data by export_figshare_fig1.m.

clear; close all
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..')); P = config_paths();
outdir = fullfile(P.out_fig1, 'chirality');
if ~exist(outdir, 'dir'); mkdir(outdir); end

all_data = load(P.framesV6).all_data;
num_samples = numel(all_data);
time_lags = get_time_lags(P, all_data);

w = 5;           % q window half-width [points]
tube_r = 4;      % centerline tube radius for the *_tube.png renders [um]
cmap = managua_br();
q_lim = [-1, 1];

%% 1. 3D renders at the first and last frame
for idx_sample = 1:num_samples
    data = all_data{idx_sample};
    x_all = arclength_resample(data.x);
    times = [1, min(size(data.x, 1), numel(data.meshes))];

    figure('Color', 'w')
    h_gut = gobjects(1,2); h_dots = gobjects(1,2); ax = gobjects(1,2);
    tube_pts = cell(1,2); tube_q = cell(1,2);
    for p = 1:2
        tidx = times(p);
        msh = data.meshes{tidx};
        curv = squeeze(x_all(tidx, 2:end, :));         % first point dropped
        [q, centres] = coiling_profile(curv, w);
        q = savgol(q, 2, 5);

        ax(p) = subplot(1, 2, p);
        h_gut(p) = trisurf(triangulation(msh.f, msh.v), 'EdgeColor', 'none', ...
            'FaceColor', 0.9*[1,1,1], 'FaceAlpha', 0.5);
        hold on
        h_dots(p) = scatter3(curv(centres,1), curv(centres,2), curv(centres,3), [], q, 'filled');
        tube_pts{p} = curv(centres,:); tube_q{p} = q(:);
        plot3([200, 250], -72*[1, 1], -60*[1, 1], 'Color', 'k', 'LineWidth', 2)   % 50 um scale bar
        colormap(cmap); clim(q_lim)
        axis equal off
        xlim([-50, 275]); ylim([-85, 85]); zlim([-75, 80])
        lighting gouraud; material dull
        view(0, -90)                                   % ventral view
        camlight; camlight(0, -90)
    end
    exportgraphics(gcf, fullfile(outdir, sprintf('chiral_sample%d.pdf', idx_sample)), ...
        'ContentType', 'image', 'Resolution', 400)

    % alternate version: centerline as a tube coloured by q, as in fig1_chiral_movie.m
    for p = 1:2
        delete(h_dots(p))
        set(h_gut(p), 'FaceAlpha', 0.35)
        [tx, ty, tz] = tubeplot(tube_pts{p}', tube_r, 16, 1e-9);   % one ring per point + end caps
        tc = repmat([tube_q{p}(1); tube_q{p}; tube_q{p}(end)]', size(tx,1), 1);
        surf(ax(p), tx, ty, tz, tc, 'EdgeColor', 'none', 'FaceColor', 'interp')
        delete(findobj(ax(p), 'Type', 'light'))
        material(ax(p), 'dull')                        % no specular glare on the tube
        camlight(ax(p)); camlight(ax(p), 0, -90)
    end
    exportgraphics(gcf, fullfile(outdir, sprintf('chiral_sample%d_tube.png', idx_sample)), 'Resolution', 600)
    close(gcf)
end

%% 2. chirality kymographs
all_measures = cell(1, num_samples);
all_s  = cell(1, num_samples);
all_tps = cell(1, num_samples);
for idx_sample = 1:num_samples
    data = all_data{idx_sample};
    folds = data.folds;
    for k = 1:3, folds(1:data.fold_onset(k), k) = -1; end

    centerlines = arclength_resample(data.x);
    num_t = size(centerlines, 1);
    % q along the centerline for every frame; 9 anterior and 1 posterior
    % point dropped; positions are s/L of the full centerline
    [measure, s_kymo] = coiling_kymograph(centerlines, w, 9, 1);

    % max-x point (anterior/posterior split) as an arclength fraction
    split_s = zeros(num_t, 1);
    for t = 1:num_t
        ce = squeeze(centerlines(t,:,:));
        ss = [0; cumsum(vecnorm(diff(ce), 2, 2))];
        [~, imax] = max(ce(:,1));
        split_s(t) = ss(imax) / ss(end);
    end

    tps = aligned_time((1:num_t)', time_lags(idx_sample), 'frame0', P);
    figure('Color', 'w')
    imagesc(s_kymo, tps, measure')
    colormap(cmap); clim(q_lim)
    hold on
    fold_s = fold_arclength(data, 'x');
    for k = 1:3
        valid = folds(:,k) > 0;
        plot(fold_s(valid,k), tps(valid), 'w-')
    end
    plot(split_s, tps, 'Color', 'w', 'LineStyle', '--')
    xlim([0.1, 0.9]); xticks(0:0.25:1)
    xlabel('s/L'); ylabel('Time [min]')
    cbar = colorbar(); cbar.Label.String = 'q';
    box on
    exportgraphics(gcf, fullfile(outdir, sprintf('kymos_%d.pdf', idx_sample)), 'ContentType', 'vector')
    close(gcf)

    all_measures{idx_sample} = measure;
    all_s{idx_sample} = s_kymo;
    all_tps{idx_sample} = tps;
end

%% 3. PCA of the kymograph rows (one observation per sample and frame)
npos_common = 100;
s_common = linspace(max(cellfun(@(s) s(1), all_s)), min(cellfun(@(s) s(end), all_s)), npos_common);
X = []; labels_sample = []; labels_time = []; labels_WT = [];
for idx_sample = 1:num_samples
    m = all_measures{idx_sample};
    nt = size(m, 2);
    m_interp = zeros(npos_common, nt);
    for tidx = 1:nt
        m_interp(:, tidx) = interp1(all_s{idx_sample}, m(:, tidx), s_common, 'linear');
    end
    X = [X; m_interp'];                                                   %#ok<AGROW>
    labels_sample = [labels_sample; idx_sample * ones(nt, 1)];            %#ok<AGROW>
    labels_time   = [labels_time; all_tps{idx_sample}(:)];                %#ok<AGROW>
    labels_WT     = [labels_WT; repmat(P.samples(idx_sample).is_WT, nt, 1)]; %#ok<AGROW>
end
labels_WT = logical(labels_WT);
X = fillmissing(X, 'constant', 0);
valid_cols = std(X, 0, 1) > eps;
[coeff_valid, score, ~, ~, explained] = pca(X(:, valid_cols));
coeff = zeros(npos_common, size(coeff_valid, 2));
coeff(valid_cols, :) = coeff_valid;

% colour by time: WT dark -> light blue, byn>Myo1C dark -> light red
t_norm = (labels_time - min(labels_time)) / (max(labels_time) - min(labels_time) + eps);
colors_time = [0.5 + 0.5*t_norm, 0.1 + 0.5*t_norm, 0.1 + 0.3*t_norm];
colors_time(labels_WT, :) = [0.1 + 0.5*t_norm(labels_WT), 0.2 + 0.5*t_norm(labels_WT), 0.5 + 0.5*t_norm(labels_WT)];
% colour by sample, as in fig1_wr.m
sample_colors = cell2mat(arrayfun(@(i) sample_color(i, num_samples), (1:num_samples)', 'UniformOutput', false));
colors_sample = sample_colors(labels_sample, :);

figure('Color', 'w')
pc_labels = {['PC1 (' num2str(explained(1), '%.1f') '%)'], ['PC2 (' num2str(explained(2), '%.1f') '%)']};
cols = {colors_time, colors_sample}; ttl = {'by time', 'by sample'};
for p = 1:2
    subplot(1, 3, p)
    scatter(score(:,1), score(:,2), 8, cols{p}, 'filled', 'MarkerFaceAlpha', 0.6)
    xlabel(pc_labels{1}); ylabel(pc_labels{2}); title(ttl{p})
    box on; axis square
end
subplot(1, 3, 3); hold on                              % legend: one swatch per sample
for i = 1:num_samples
    fill([0, 1, 1, 0], [i-0.5, i-0.5, i+0.5, i+0.5], sample_colors(i, :), 'EdgeColor', 'none')
    if P.samples(i).is_WT
        lbl = sprintf('WT %d', sum([P.samples(1:i).is_WT]));
    else
        lbl = sprintf('Myo %d', sum(~[P.samples(1:i).is_WT]));
    end
    text(1.15, i, lbl, 'FontSize', 7, 'VerticalAlignment', 'middle')
end
xlim([0, 3]); ylim([0.5, num_samples + 0.5])
set(gca, 'YDir', 'reverse'); axis off
exportgraphics(gcf, fullfile(outdir, 'pca_scatter.pdf'), 'ContentType', 'vector')

% PC1 and PC2 loadings as 1-D images
figure('Color', 'w')
for pc = 1:2
    subplot(2, 1, pc)
    imagesc(s_common, 1, coeff(:,pc)')
    colormap(cmap)
    clim(max(abs(coeff(:,pc))) * [-1, 1])
    colorbar()
    xlabel('s/L'); ylabel(['PC' num2str(pc)]); yticks([])
    box on
end
exportgraphics(gcf, fullfile(outdir, 'pca_loadings.pdf'), 'ContentType', 'vector')
