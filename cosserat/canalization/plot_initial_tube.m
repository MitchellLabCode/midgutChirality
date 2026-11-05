
% load reference mesh
path_meshes = "../data/meshes_spRS/meshes_sample6.mat";
path_data = "../data/refIV6_sample6_walls3_cball.mat";


path_meshes = "../data/meshes_spRS/meshes_sample1.mat";
path_data = "../data/refIV6_sample1_walls3_cball.mat";
%path_meshes = "../data/meshes_spRS/meshes_sample4.mat";
%path_data = "../data/refIV6_sample4_walls3_cball.mat";

meshes = load(path_meshes);
meshes = meshes.meshes;
data = load(path_data);

data_x = data.x;
data_Q = data.Q;
num_t = size(data_x, 1);




% load simulation ouputs

%folder_name = "20260209_Q_Tp_t45";
%folder_name = "20260218_Q_twp_t54_s6V";
t0 = 54;
folder_name = "20260325_canalization_box2f_s6s3"; % with fixed posterior IC

folder_name = "20260406_canalization_box2f_nocontact_notail_s6s1"; %

folder_name = "20260402_canalization_box2f_s6s1"; %


%path_sim = fullfile(folder_name, "cutoff01_T1.0_p0.3.mat");
%sim_data = load(path_sim);


% get maximal box
xmin = 1e9;
xmax = -1e9;
ymax = -1e9;
ymin = 1e9;
zmin = 1e9;
zmax = -1e9;
for t = t0:num_t
    verts = meshes{t}.v;
    mx = min(verts(:,1));
    Mx = max(verts(:,1));
    my = min(verts(:,2));
    My = max(verts(:,2));
    mz= min(verts(:,3));
    Mz = max(verts(:,3));

    if mx < xmin
        xmin = mx;
    end
    if Mx > xmax
        xmax = Mx;
    end
    if my < ymin
        ymin = my;
    end
    if My > ymax
        ymax = My;
    end
    if mz < zmin
        zmin = mz;
    end
    if Mz > zmax
        zmax = Mz;
    end
end


vik = load("vik.mat").vik;

verts = meshes{num_t}.v;
faces = meshes{num_t}.f;

% interpolate data to get centerline points at the same time steps as the meshes
%xsaves = sim_data.xsaves;
%num_tsim = size(xsaves, 1);
%x_interp = interp1(linspace(t0, num_t, num_tsim), squeeze(xsaves(:,:,1)), linspace(t0, num_t, num_t-t0+1));
%y_interp = interp1(linspace(t0, num_t, num_tsim), squeeze(xsaves(:,:,2)), linspace(t0, num_t, num_t-t0+1));
%z_interp = interp1(linspace(t0, num_t, num_tsim), squeeze(xsaves(:,:,3)), linspace(t0, num_t, num_t-t0+1));

% rescale x, y, z by length of the rod at t0
L_t0 = sum(sqrt(sum(diff(squeeze(data_x(t0,:,:)), 1, 1).^2, 2)));
%x_interp = x_interp * L_t0;
%y_interp = y_interp * L_t0;
%z_interp = z_interp * L_t0;

figure
mat_files = dir(fullfile(folder_name, '*.mat'));
disp(['Found ' num2str(length(mat_files)) ' .mat files in folder: ' folder_name]);
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
num_files = length(mat_files)/2;
figure;
num_cols = 6; % Number of columns for subplots
num_rows = ceil(2*num_files / num_cols); % Calculate number of rows needed

x0=10;
y0=10;
width=1000;
height=1000;
set(gcf,'position',[x0,y0,width,height])

tiledlayout(num_rows, num_cols, 'Padding', 'tight', 'TileSpacing', 'compact'); 
colormap(vik)

for k = 1:2:length(mat_files)
    % Load each .mat file
    sim_data = load(fullfile(folder_name, mat_files(k).name));
    n = n_values(k);
    tau = tau_values(k);
    xsaves = sim_data.xsaves;
    Qsaves = sim_data.Qsaves;

    d1 = squeeze(data_Q(1, :, 1,:));
    d1_sim = squeeze(Qsaves(1, :, 1,:));

    %size()
    %subplot(num_rows, num_cols, k);
    nexttile;
    %trisurf(triangulation(faces, verts), 'edgecolor', 'none', ...
    %        'facecolor', 0.9* [1, 1, 1], 'faceAlpha', 0.5);
    hold on
    % plot centerline ref + sim
    %plot3(data_x(num_t,:,1), data_x(num_t,:,2), data_x(num_t,:,3), 'k-', 'linewidth', 2)
    %plot3(xsaves(end,:, 1)* L_t0, xsaves(end,:, 2)* L_t0, xsaves(end,:, 3)* L_t0, 'r-', 'linewidth', 2)

    sim_x = [xsaves(1,:, 1)* L_t0; xsaves(1,:, 2)* L_t0; xsaves(1,:, 3)* L_t0];
    %disp(size(sim_x))
    [xx, yy, zz] = tubeplot(sim_x,  15, 24);
    %surf(xx,yy,zz, 'Facecolor', [212, 99, 74]/255, 'EdgeColor','none')

    surf(xx,yy,zz, tau *ones(size(zz)), 'EdgeColor','none')
    clim([min(tau_values), max(tau_values)])


    %size(xsaves(end,:))

    %quiver3(squeeze(data_x(num_t,1:end-1,1))', squeeze(data_x(num_t,1:end-1,2))', squeeze(data_x(num_t,1:end-1,3))', squeeze(d1(:,1)), squeeze(d1(:,2)), squeeze(d1(:,3)),  'k-', 'linewidth', 2, 'AutoScale', 'on', 'AutoScaleFactor', 0.5)
    %quiver3(squeeze(xsaves(end,1:end-1, 1)* L_t0)', squeeze(xsaves(end,1:end-1, 2)* L_t0)', squeeze(xsaves(end,1:end-1, 3)* L_t0)', squeeze(d1_sim(:,1)), squeeze(d1_sim(:,2)), squeeze(d1_sim(:,3)),'r-', 'linewidth', 2, 'AutoScale', 'on', 'AutoScaleFactor', 0.5)

    %axis([xmin, xmax, ymin, ymax, zmin, zmax])
    daspect([1,1,1])
    axis off

    xlim tight
    ylim tight
    zlim tight

    title(['n = ' num2str(n) ', tau = ' num2str(tau)], 'fontsize', 6)

    lighting gouraud
    material dull
    view(0,-90) % ventral view   
    %camlight(0,90)
    camlight("headlight")

    %view(0,-90) % ventral view    


    nexttile;
    %trisurf(triangulation(faces, verts), 'edgecolor', 'none', ...
    %        'facecolor', 0.9* [1, 1, 1], 'faceAlpha', 0.5);
    hold on
    % plot centerline ref + sim
    %plot3(data_x(num_t,:,1), data_x(num_t,:,2), data_x(num_t,:,3), 'k-', 'linewidth', 2)
    %plot3(xsaves(end,:, 1)* L_t0, xsaves(end,:, 2)* L_t0, xsaves(end,:, 3)* L_t0, 'r-', 'linewidth', 2)

    sim_x = [xsaves(1,:, 1)* L_t0; xsaves(1,:, 2)* L_t0; xsaves(1,:, 3)* L_t0];
    %disp(size(sim_x))
    [xx, yy, zz] = tubeplot(sim_x,  15, 24);
    %surf(xx,yy,zz, 'Facecolor', [212, 99, 74]/255, 'EdgeColor','none')
    surf(xx,yy,zz, tau *ones(size(zz)), 'EdgeColor','none')
    clim([min(tau_values), max(tau_values)])

    % quiver3(squeeze(data_x(num_t,1:end-1,1))', squeeze(data_x(num_t,1:end-1,2))', squeeze(data_x(num_t,1:end-1,3))', squeeze(d1(:,1)), squeeze(d1(:,2)), squeeze(d1(:,3)),  'k-', 'linewidth', 2, 'AutoScale', 'on', 'AutoScaleFactor', 0.5)
    % quiver3(squeeze(xsaves(end,1:end-1, 1)* L_t0)', squeeze(xsaves(end,1:end-1, 2)* L_t0)', squeeze(xsaves(end,1:end-1, 3)* L_t0)', squeeze(d1_sim(:,1)), squeeze(d1_sim(:,2)), squeeze(d1_sim(:,3)),'r-', 'linewidth', 2, 'AutoScale', 'on', 'AutoScaleFactor', 0.5)


    %axis([xmin, xmax, ymin, ymax, zmin, zmax])
    daspect([1,1,1])

    xlim tight
    ylim tight
    zlim tight

    lighting gouraud
    material dull

    
    axis off
    %title(['T = ' num2str(T) ', p = ' num2str(p)], 'fontsize', 6)

    camlight(0,-90)
    view([0,0]) % side view 

    %view([60,40])% perspective view
    %view(45,30) % perspective view
    %saveas(gcf, fullfile(folder_name, "images/perspective/", ['frame' num2str(t) '.png']))
    %close(gcf)
end
save_filename = [char(folder_name) '_twin_tube_half_initial.png'];
%saveas(gcf, save_filename)
print(save_filename, '-dpng', '-r500') % Save at 500 DPI

    
