% Overlay version of plot_tube_canal.m: final shapes of 4 runs, overlaid in
% pairs with transparency. Columns: the two alpha < 0 runs, the two alpha > 0
% runs. Rows: ventral and side views, on common axes so the panels compare.
% alpha = 1 - 2*tau, as in plot_canalization_frames.m.
% Colours as in plot_tube_canal.m: tau on vik over [min tau, max tau]
% (= alpha on vik_r over [-1.2, 1.2]).

% load reference data (for the length scale L_t0 only)
path_data = "../data/refIV6_sample1_walls3_cball.mat";
%path_data = "../data/refIV6_sample6_walls3_cball.mat";
data = load(path_data);
data_x = data.x;
t0 = 54;

folder_name = "202609022_canalization_box2f_s6s1";

alpha_pairs = {[-1, -0.5], [0.5, 1]};   % one panel per pair (snapped to the nearest run)
face_alpha  = 0.55;                      % transparency of the overlaid tubes
tube_r      = 15;                        % as in plot_tube_canal.m (data units)

% rescale x, y, z by length of the rod at t0
L_t0 = sum(sqrt(sum(diff(squeeze(data_x(t0,:,:)), 1, 1).^2, 2)));

vik = load("vik.mat").vik;

mat_files = dir(fullfile(folder_name, '*.mat'));
disp(['Found ' num2str(length(mat_files)) ' .mat files in folder: ' char(folder_name)]);
% Sort mat_files by n and tau values
n_values = nan(length(mat_files), 1);
tau_values = nan(length(mat_files), 1);

for k = 1:length(mat_files)
    filename = mat_files(k).name;
    tokens = regexp(filename, 'n([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)_tau([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)\.mat', 'tokens');
    if ~isempty(tokens)
        n_values(k) = str2double(tokens{1}{1});
        tau_values(k) = str2double(tokens{1}{2});
    end
end

valid_idx = ~isnan(n_values) & ~isnan(tau_values);
mat_files = mat_files(valid_idx);
n_values = n_values(valid_idx);
tau_values = tau_values(valid_idx);

[~, sort_idx] = sortrows([n_values, tau_values]);
mat_files = mat_files(sort_idx);
n_values = n_values(sort_idx);
tau_values = tau_values(sort_idx);
alpha_values = 1 - 2*tau_values;

% colour of a run, as given by the surf CData = tau and clim below (for the titles)
tmin = min(tau_values); tmax = max(tau_values); nv = size(vik, 1);
col_of = @(tau) vik(1 + round((tau - tmin)/(tmax - tmin)*(nv-1)), :);

% select the runs nearest to the requested alpha
sel = cell(size(alpha_pairs));
for c = 1:numel(alpha_pairs)
    for a = alpha_pairs{c}
        [~, m] = min(abs(alpha_values - a));
        fprintf('alpha %g -> %.3f (%s)\n', a, alpha_values(m), mat_files(m).name);
        sel{c}(end+1) = m;
    end
end

% final shapes, and one box over all of them
shapes = cell(length(mat_files), 1);
b = [inf -inf inf -inf inf -inf];
for m = [sel{:}]
    sim_data = load(fullfile(folder_name, mat_files(m).name));
    xsaves = sim_data.xsaves;
    sim_x = [xsaves(end,:, 1)* L_t0; xsaves(end,:, 2)* L_t0; xsaves(end,:, 3)* L_t0];
    shapes{m} = sim_x;
    b = [min(b(1), min(sim_x(1,:))) max(b(2), max(sim_x(1,:))) ...
         min(b(3), min(sim_x(2,:))) max(b(4), max(sim_x(2,:))) ...
         min(b(5), min(sim_x(3,:))) max(b(6), max(sim_x(3,:)))];
end
box = b + 1.5*tube_r*[-1 1 -1 1 -1 1];

figure;
x0=10;
y0=10;
width=1000;
height=700;
set(gcf,'position',[x0,y0,width,height])

views = {'ventral', 'side'};
tiledlayout(numel(views), numel(alpha_pairs), 'Padding', 'tight', 'TileSpacing', 'compact');
colormap(vik)

for iv = 1:numel(views)
    for c = 1:numel(alpha_pairs)
        nexttile;
        hold on
        lab = {};
        for m = sel{c}
            [xx, yy, zz] = tubeplot(shapes{m}, tube_r, 24);
            surf(xx,yy,zz, tau_values(m) *ones(size(zz)), 'EdgeColor','none', 'FaceAlpha', face_alpha)
            lab{end+1} = sprintf('\\color[rgb]{%.3f,%.3f,%.3f}\\alpha = %.2f', col_of(tau_values(m)), alpha_values(m));
        end
        clim([min(tau_values), max(tau_values)])

        axis(box)
        daspect([1,1,1])
        axis off

        lighting gouraud
        material dull
        if strcmp(views{iv}, 'side')
            camlight(0,-90)
            view([0,0]) % side view
        else
            view(0,-90) % ventral view
            camlight("headlight")
        end

        if iv == 1
            title({views{iv}, strjoin(lab, '    ')}, 'fontsize', 9)
        else
            title(views{iv}, 'fontsize', 9)
        end
    end
end
save_filename = [char(folder_name) '_twin_tube_overlay.png'];
print(save_filename, '-dpng', '-r500') % Save at 500 DPI
