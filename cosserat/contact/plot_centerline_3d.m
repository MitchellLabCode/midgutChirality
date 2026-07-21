
% load reference mesh
path_meshes = "../data/meshes_spRS/meshes_sample3.mat";
path_data = "../data/refIV5_sample3_walls3.mat";

meshes = load(path_meshes);
meshes = meshes.meshes;
data = load(path_data);

data_x = data.x;
num_t = size(data_x, 1);

t0 = 50;

t0 = 54;

t0 = 48;
% load simulation ouputs

folder_name = "./20260216_Q_Tp_finer_t54_s6V";
folder_name = "./20260220_Q_twp_t54f2_s6V";

folder_name = "./20260304_Q_twp_contact_t48_w3_a1e2_s3V5";

path_sim = fullfile(folder_name, "cutoff01_T3.0_p10.0.mat");
path_sim = fullfile(folder_name, "cutoff01contact_tw3.0_p10.0.mat");
sim_data = load(path_sim);


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
    mz = min(verts(:,3));
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

% print all frames in a loop


if ~isfolder(fullfile(folder_name, "images"))
    mkdir(fullfile(folder_name, "images"))
end
if ~isfolder(fullfile(folder_name, "images/dorsal"))
    mkdir(fullfile(folder_name, "images/dorsal"))
end
if ~isfolder(fullfile(folder_name, "images/perspective"))
    mkdir(fullfile(folder_name, "images/perspective"))
end
if ~isfolder(fullfile(folder_name, "images/side"))
    mkdir(fullfile(folder_name, "images/side"))
end

% interpolate data to get centerline points at the same time steps as the meshes
xsaves = sim_data.xsaves;
num_tsim = size(xsaves, 1);
x_interp = interp1(linspace(t0, num_t, num_tsim), squeeze(xsaves(:,:,1)), linspace(t0, num_t, num_t-t0+1));
y_interp = interp1(linspace(t0, num_t, num_tsim), squeeze(xsaves(:,:,2)), linspace(t0, num_t, num_t-t0+1));
z_interp = interp1(linspace(t0, num_t, num_tsim), squeeze(xsaves(:,:,3)), linspace(t0, num_t, num_t-t0+1));

% rescale x, y, z by length of the rod at t0
L_t0 = sum(sqrt(sum(diff(squeeze(data_x(t0,:,:)), 1, 1).^2, 2)));
x_interp = x_interp * L_t0;
y_interp = y_interp * L_t0;
z_interp = z_interp * L_t0;

figure
for t = t0:num_t
    verts = meshes{t}.v;
    faces = meshes{t}.f;
    clf;
    trisurf(triangulation(faces, verts), 'edgecolor', 'none', ...
            'facecolor', 0.9* [1, 1, 1], 'faceAlpha', 0.5);
    hold on
    % plot centerline ref + sim
    plot3(data_x(t,:,1), data_x(t,:,2), data_x(t,:,3), 'k-', 'linewidth', 2)
    plot3(x_interp(t-t0+1,:), y_interp(t-t0+1,:), z_interp(t-t0+1,:), 'r-', 'linewidth', 2)

    axis([xmin, xmax, ymin, ymax, zmin, zmax])
    daspect([1,1,1])
    axis off

    view(0,90) % dorsal view
    saveas(gcf, fullfile(folder_name, "images/dorsal/", ['frame' num2str(t) '.png']))

    view([60,40])% perspective view
    %view(45,30) % perspective view
    saveas(gcf, fullfile(folder_name, "images/perspective/", ['frame' num2str(t) '.png']))

    view([0,0])% side view (from left)
    saveas(gcf, fullfile(folder_name, "images/side/", ['frame' num2str(t) '.png']))
    %close(gcf)
end

    
