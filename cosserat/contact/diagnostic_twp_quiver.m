
% load reference mesh


%path_meshes = "../data/meshes_spRS/meshes_sample1.mat";
%path_data = "../data/refIV5_sample1_walls3_cball.mat";

%path_meshes = "../data/meshes_spRS/meshes_sample2.mat";
%path_data = "../data/refIV5_sample2_walls.mat";

%path_meshes = "../data/meshes_spRS/meshes_sample5.mat";
%path_data = "../data/refIV5_sample5_walls3_cball.mat";

%path_meshes = "../data/meshes_spRS/meshes_sample4.mat";
%path_data = "../data/ref_sample4_walls.mat";


path_meshes = "../data/meshes_spRS/meshes_sample7.mat";
path_data = "../data/refIV5_sample7_walls3_cball.mat";

meshes = load(path_meshes);
meshes = meshes.meshes;
data = load(path_data);

data_x = data.x;
data_Q = data.Q;
num_t = size(data_x, 1);




% load simulation ouputs

t0 = 45;
folder_name = "20260304_Q_twp_contactf_fact01_t48_w3_a5e2_R1d2_s3V5"; 
folder_name = "20260304_Q_twp_contact_t45f_w3_a3e2_R1d0_s3V5";

t0 = 36;
folder_name = "20260305_Q_twp_contactf_fact01_t48_w3_a5e2_R1d2_s3V5"; 
folder_name = "20260311_Q_twp_contact_t36_w3_a3e2_R1d0_s4V5";

t0 = 32;
folder_name = '20260311_Q_twp_contact_t32f_w3_a3e2_R1d0_s4V5';

t0 =25;
folder_name = '20260506_Q_twp_contact_t25_w3_a3e2_R1d0_s2V5';


%t0 = 38;
%folder_name = "20260402_Q_twp_contact_t38a0_w3_a3e2_R1d0_s1V5";

%t0 = 25;
%folder_name = "20260505_Q_twp_contact_t25_w3_a3e2_R1d0_s5V5";


%t0 = 54;
%folder_name = "20260304_Q_twp_contact_t54f_w3_a1e2_nl_s6V5"; 

t0 =25;
folder_name = "20260708_Q_twp_contact_t25_w3_a3e2_R1d0_s7V5"; 



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
% Sort mat_files by T and p values
T_values = [];
p_values = [];

for k = 1:length(mat_files)
    filename = mat_files(k).name;
    tokens = regexp(filename, 'tw([\d.]+)_p([\d.]+)', 'tokens');
    if ~isempty(tokens)
        T_values(k) = str2double(tokens{1}{1});
        p_values(k) = str2double(tokens{1}{2}(1:end-1));
    end
end


[~, sort_idx] = sortrows([T_values', p_values']);
mat_files = mat_files(sort_idx);
num_files = length(mat_files);
figure;
num_cols = 6; % Number of columns for subplots
num_rows = ceil(2*num_files / num_cols); % Calculate number of rows needed

x0=10;
y0=10;
width=1000;
height=1000;
set(gcf,'position',[x0,y0,width,height])

tiledlayout(num_rows, num_cols, 'Padding', 'tight', 'TileSpacing', 'compact'); 

for k = 1:length(mat_files)
    % Load each .mat file
    sim_data = load(fullfile(folder_name, mat_files(k).name));
    % Extract T and p from filename
    filename = mat_files(k).name;

    tokens = regexp(filename, 'tw([\d.]+)_p([\d.]+)', 'tokens');
    if ~isempty(tokens)
        T = str2double(tokens{1}{1});
        p =tokens{1}{2};
        p = str2double(p(1:end-1));
    end
    xsaves = sim_data.xsaves;
    Qsaves = sim_data.Qsaves;

    d1 = squeeze(data_Q(end, :, 1,:));
    d1_sim = squeeze(Qsaves(end, :, 1,:));

    %size()
    %subplot(num_rows, num_cols, k);
    nexttile;
    trisurf(triangulation(faces, verts), 'edgecolor', 'none', ...
            'facecolor', 0.9* [1, 1, 1], 'faceAlpha', 0.5);
    hold on
    % plot centerline ref + sim
    plot3(data_x(num_t,:,1), data_x(num_t,:,2), data_x(num_t,:,3), 'k-', 'linewidth', 2)
    plot3(xsaves(end,:, 1)* L_t0, xsaves(end,:, 2)* L_t0, xsaves(end,:, 3)* L_t0, 'r-', 'linewidth', 2)

    %size(xsaves(end,:))

    %quiver3(squeeze(data_x(num_t,1:end-1,1))', squeeze(data_x(num_t,1:end-1,2))', squeeze(data_x(num_t,1:end-1,3))', squeeze(d1(:,1)), squeeze(d1(:,2)), squeeze(d1(:,3)),  'k-', 'linewidth', 2, 'AutoScale', 'on', 'AutoScaleFactor', 0.5)
    %quiver3(squeeze(xsaves(end,1:end-1, 1)* L_t0)', squeeze(xsaves(end,1:end-1, 2)* L_t0)', squeeze(xsaves(end,1:end-1, 3)* L_t0)', squeeze(d1_sim(:,1)), squeeze(d1_sim(:,2)), squeeze(d1_sim(:,3)),'r-', 'linewidth', 2, 'AutoScale', 'on', 'AutoScaleFactor', 0.5)

    axis([xmin, xmax, ymin, ymax, zmin, zmax])
    daspect([1,1,1])
    axis off
    title(['tw = ' num2str(T) ', p = ' num2str(p)], 'fontsize', 6)

    view(0,90) % dorsal view    


    nexttile;
    trisurf(triangulation(faces, verts), 'edgecolor', 'none', ...
            'facecolor', 0.9* [1, 1, 1], 'faceAlpha', 0.5);
    hold on
    % plot centerline ref + sim
    plot3(data_x(num_t,:,1), data_x(num_t,:,2), data_x(num_t,:,3), 'k-', 'linewidth', 2)
    plot3(xsaves(end,:, 1)* L_t0, xsaves(end,:, 2)* L_t0, xsaves(end,:, 3)* L_t0, 'r-', 'linewidth', 2)

    %quiver3(squeeze(data_x(num_t,1:end-1,1))', squeeze(data_x(num_t,1:end-1,2))', squeeze(data_x(num_t,1:end-1,3))', squeeze(d1(:,1)), squeeze(d1(:,2)), squeeze(d1(:,3)),  'k-', 'linewidth', 2, 'AutoScale', 'on', 'AutoScaleFactor', 0.5)
    %quiver3(squeeze(xsaves(end,1:end-1, 1)* L_t0)', squeeze(xsaves(end,1:end-1, 2)* L_t0)', squeeze(xsaves(end,1:end-1, 3)* L_t0)', squeeze(d1_sim(:,1)), squeeze(d1_sim(:,2)), squeeze(d1_sim(:,3)),'r-', 'linewidth', 2, 'AutoScale', 'on', 'AutoScaleFactor', 0.5)


    axis([xmin, xmax, ymin, ymax, zmin, zmax])
    daspect([1,1,1])
    axis off
    %title(['T = ' num2str(T) ', p = ' num2str(p)], 'fontsize', 6)

    view([0,0]) % side view 

    %view([60,40])% perspective view
    %view(45,30) % perspective view
    %saveas(gcf, fullfile(folder_name, "images/perspective/", ['frame' num2str(t) '.png']))
    %close(gcf)
end
save_filename = [char(folder_name) '_twin_noquiver.png'];
%saveas(gcf, save_filename)
print(save_filename, '-dpng', '-r500') % Save at 500 DPI

    
