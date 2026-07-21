
% load reference mesh
path_meshes = "../data/meshes_spRS/meshes_sample2.mat";
path_data = "../data/refIV5_sample2_walls3_cball.mat";

addpath(".")
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


t0 = 37;
folder_name = "./20260331_Q_twp_contact_t37_w3_a3e2_R1d0_s1V5";

t0 = 30;
folder_name = "./20260310_Q_twp_contact_t30fa2_w3_a3e2_R1d0_s4V5";


t0 = 25;
folder_name = "./20260505_Q_twp_contact_t25_w3_a3e2_R1d0_s5V5";

t0 = 32;
folder_name = '20260311_Q_twp_contact_t32f_w3_a3e2_R1d0_s4V5'

t0 =25;
folder_name = '20260506_Q_twp_contact_t25_w3_a3e2_R1d0_s2V5'

path_sim = fullfile(folder_name, "cutoff01contact_tw3.0_p10.0.mat");
sim_data = load(path_sim);


vik = load('vik.mat').vik;

vik_b = vik(22,:); % WT
vik_r = vik(234,:); % Myo

vik_c = vik_b; % set which color to use


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


if ~isfolder(fullfile(folder_name, "imagesQtube"))
    mkdir(fullfile(folder_name, "imagesQtube"))
end
if ~isfolder(fullfile(folder_name, "imagesQtube/dorsal"))
    mkdir(fullfile(folder_name, "imagesQtube/dorsal"))
end
if ~isfolder(fullfile(folder_name, "imagesQtube/ventral"))
    mkdir(fullfile(folder_name, "imagesQtube/ventral"))
end
if ~isfolder(fullfile(folder_name, "imagesQtube/perspective"))
    mkdir(fullfile(folder_name, "imagesQtube/perspective"))
end
if ~isfolder(fullfile(folder_name, "imagesQtube/side"))
    mkdir(fullfile(folder_name, "imagesQtube/side"))
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
            'facecolor', 0.9* [1, 1, 1], 'faceAlpha', 0.2);
    hold on
    % plot centerline ref + sim
    %[xx, yy, zz] = tubeplot(squeeze(data_x(t,:,:))',  1.5, 12);
    [xx, yy, zz] = tubeplot(squeeze(data_x(t,:,:))',  5, 24);
    surf(xx,yy,zz, 'Facecolor', 0.6* [1, 1, 1], 'EdgeColor','none')
    %plot3(data_x(t,:,1), data_x(t,:,2), data_x(t,:,3), 'k-', 'linewidth', 2)
    sim_x = [squeeze(x_interp(t-t0+1,:)); squeeze(y_interp(t-t0+1,:)); squeeze(z_interp(t-t0+1,:))];
    %disp(size(sim_x))
    %[xx, yy, zz] = tubeplot(sim_x,  1.5, 12);
    [xx, yy, zz] = tubeplot(sim_x,  5, 24);
    %surf(xx,yy,zz, 'Facecolor', [212, 99, 74]/255, 'EdgeColor','none')
    surf(xx,yy,zz, 'Facecolor', vik_c, 'EdgeColor','none')
    %plot3(x_interp(t-t0+1,:), y_interp(t-t0+1,:), z_interp(t-t0+1,:), 'r-', 'linewidth', 2)

    %size(data_x(t,1:end-1,1))
    %size(d1(:,1))

    quiver3(squeeze(data_x(t,1:end-1,1))', squeeze(data_x(t,1:end-1,2))', squeeze(data_x(t,1:end-1,3))', squeeze(d1(:,1)), squeeze(d1(:,2)), squeeze(d1(:,3)), 'k-', 'linewidth', 2, 'AutoScale', 'on', 'AutoScaleFactor', 0.4)

    %size(x_interp)
    %size(d1x_interp)
    %quiver3(x_interp(t-t0+1,1:end-1), y_interp(t-t0+1,1:end-1), z_interp(t-t0+1,1:end-1), d1x_interp(t-t0+1,:), d1y_interp(t-t0+1,:), d1z_interp(t-t0+1,:), 'r-', 'linewidth', 2, 'AutoScale', 'on', 'AutoScaleFactor', 0.4)
    quiver3(x_interp(t-t0+1,1:end-1), y_interp(t-t0+1,1:end-1), z_interp(t-t0+1,1:end-1), d1x_interp(t-t0+1,:), d1y_interp(t-t0+1,:), d1z_interp(t-t0+1,:), 'color', vik_c, 'linewidth', 2, 'AutoScale', 'on', 'AutoScaleFactor', 0.4)

    axis([xmin, xmax, ymin, ymax, zmin, zmax])
    daspect([1,1,1])
    lighting gouraud
    material dull
    %camlight(0, -90)
    axis off

    view(0,90) % dorsal view
    cl = camlight('headlight');
    saveas(gcf, fullfile(folder_name, "imagesQtube/dorsal/", ['frame' num2str(t) '.png']))

    view(0,-90) % ventral view
    camlight(cl, 'headlight');
    saveas(gcf, fullfile(folder_name, "imagesQtube/ventral/", ['frame' num2str(t) '.png']))

    view([60,40])% perspective view
    %view(45,30) % perspective view
    camlight(cl, 'headlight');
    saveas(gcf, fullfile(folder_name, "imagesQtube/perspective/", ['frame' num2str(t) '.png']))
    %close(gcf)

    view([0,0])% side view (from left)
    camlight(cl, 'headlight');
    saveas(gcf, fullfile(folder_name, "imagesQtube/side/", ['frame' num2str(t) '.png']))
end

    
