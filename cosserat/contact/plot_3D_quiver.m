
% load reference mesh
path_meshes = "../data/meshes_spRS/meshes_sample6.mat";
path_data = "../data/refIV5_sample6_walls3.mat";

meshes = load(path_meshes);
meshes = meshes.meshes;
data = load(path_data);

data_x = data.x;
data_Q = data.Q;
num_t = size(data_x, 1);

t0 = 45;


% load simulation ouputs

folder_name = "./20260304_Q_twp_contact_t45f_w3_a3e2_R1d0_s3V5";

t0 = 54;

folder_name = "./20260304_Q_twp_contact_t54f_w3_a1e2_nl_s6V5";

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


if ~isfolder(fullfile(folder_name, "imagesQ"))
    mkdir(fullfile(folder_name, "imagesQ"))
end
if ~isfolder(fullfile(folder_name, "imagesQ/dorsal"))
    mkdir(fullfile(folder_name, "imagesQ/dorsal"))
end
if ~isfolder(fullfile(folder_name, "imagesQ/perspective"))
    mkdir(fullfile(folder_name, "imagesQ/perspective"))
end
if ~isfolder(fullfile(folder_name, "imagesQ/side"))
    mkdir(fullfile(folder_name, "imagesQ/side"))
end

% interpolate data to get centerline points at the same time steps as the meshes
xsaves = sim_data.xsaves;
Qsaves = sim_data.Qsaves;
num_tsim = size(xsaves, 1);
x_interp = interp1(linspace(t0, num_t, num_tsim), squeeze(xsaves(:,:,1)), linspace(t0, num_t, num_t-t0+1));
y_interp = interp1(linspace(t0, num_t, num_tsim), squeeze(xsaves(:,:,2)), linspace(t0, num_t, num_t-t0+1));
z_interp = interp1(linspace(t0, num_t, num_tsim), squeeze(xsaves(:,:,3)), linspace(t0, num_t, num_t-t0+1));

d1x_interp = interp1(linspace(t0, num_t, num_tsim), squeeze(Qsaves(:,:,1,1)), linspace(t0, num_t, num_t-t0+1));
d1y_interp = interp1(linspace(t0, num_t, num_tsim), squeeze(Qsaves(:,:,1,2)), linspace(t0, num_t, num_t-t0+1));
d1z_interp = interp1(linspace(t0, num_t, num_tsim), squeeze(Qsaves(:,:,1,3)), linspace(t0, num_t, num_t-t0+1));


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

    d1 = squeeze(data_Q(t,:,1,:));
    trisurf(triangulation(faces, verts), 'edgecolor', 'none', ...
            'facecolor', 0.9* [1, 1, 1], 'faceAlpha', 0.3);
    hold on
    % plot centerline ref + sim
    plot3(data_x(t,:,1), data_x(t,:,2), data_x(t,:,3), 'k-', 'linewidth', 2)
    plot3(x_interp(t-t0+1,:), y_interp(t-t0+1,:), z_interp(t-t0+1,:), 'r-', 'linewidth', 2)

    %size(data_x(t,1:end-1,1))
    %size(d1(:,1))

    quiver3(squeeze(data_x(t,1:end-1,1))', squeeze(data_x(t,1:end-1,2))', squeeze(data_x(t,1:end-1,3))', squeeze(d1(:,1)), squeeze(d1(:,2)), squeeze(d1(:,3)), 'k-', 'linewidth', 2, 'AutoScale', 'on', 'AutoScaleFactor', 0.3)

    %size(x_interp)
    %size(d1x_interp)
    quiver3(x_interp(t-t0+1,1:end-1), y_interp(t-t0+1,1:end-1), z_interp(t-t0+1,1:end-1), d1x_interp(t-t0+1,:), d1y_interp(t-t0+1,:), d1z_interp(t-t0+1,:), 'r-', 'linewidth', 2, 'AutoScale', 'on', 'AutoScaleFactor', 0.3)

    axis([xmin, xmax, ymin, ymax, zmin, zmax])
    daspect([1,1,1])
    material dull
    lighting gouraud
    camlight %(0, -90)
    axis off

    view(0,90) % dorsal view
    saveas(gcf, fullfile(folder_name, "imagesQ/dorsal/", ['frame' num2str(t) '.png']))

    view([60,40])% perspective view
    %view(45,30) % perspective view
    saveas(gcf, fullfile(folder_name, "imagesQ/perspective/", ['frame' num2str(t) '.png']))
    %close(gcf)

    view([0,0])% side view (from left)
    saveas(gcf, fullfile(folder_name, "imagesQ/side/", ['frame' num2str(t) '.png']))
end

    
