% Frame-by-frame tube renderings of every run in selected/, written as PNGs.
%
% For each selected/*.mat, renders one PNG per reference frame t0..num_t with
% the same settings as plot_tube_dynamics.m: ventral and side views side by
% side, embryo mesh (faint), imaged rod (grey), simulated rod (sample colour),
% and d1 material frames (black = imaged, colour = simulation). The images are
% numbered sequentially (frame_001.png, ...) so they stack straight into a
% movie (ImageJ, ffmpeg -i frame_%03d.png ...).
%
% Output: selected/frames/<run name>/frame_NNN.png
%
% Run from prod/. No local functions, so it can be pasted into the
% command window.

addpath(".")

% ======================= CONFIG ============================================
selected_dir = 'selected';
out_root     = fullfile(selected_dir, 'frames');

% The selected files do not record t0 (it is in the name of the source folder
% <date>_QTpTw_contact2_t<t0>_w3_a3e2_s<N>V6). These are the t0 used by
% make_selected.sh; s7 and s10 start from t0 = 20 and 35.
t0_of = containers.Map('KeyType','double','ValueType','double');
t0_of(1) = 38;  t0_of(2) = 25;  t0_of(3) = 45;  t0_of(4) = 36;  t0_of(5) = 25;
t0_of(6) = 54;  t0_of(7) = 20;  t0_of(8) = 37;  t0_of(9) = 45;  t0_of(10) = 35;

dt_min = 2;                      % minutes between reference frames

% true : one bounding box shared by ALL selected samples (frames comparable
%        across samples). false : one fixed box per sample. Either way the box
%        is fixed within a sample's sequence, so the frames do not jitter.
shared_across_samples = true;
fit_simulation = true;           % include the simulated centerline in the box
axes_pad = 10;                   % microns

tube_r = 5;  tube_n = 24;
fig_px = [1400 650];             % fixed pixel size -> every PNG identical size
save_dpi = '-r150';

draw_quivers   = true;
director_idx   = 1;
quiver_len     = 25;             % microns (fixed; AutoScale off)
quiver_stride  = 6;
quiver_lw      = 1.0;
draw_ref_frame = true;
% ===========================================================================

vik = load('vik.mat').vik;

files = dir(fullfile(selected_dir, 'sample*_*_T*_tw*_p*.mat'));
if isempty(files)
    error('No runs found in %s', selected_dir);
end

% ---- parse names; load everything needed for the bounding boxes -----------
runs = struct('file',{},'sid',{},'geno',{},'T',{},'tw',{},'p',{},'t0',{});
for k = 1:numel(files)
    tok = regexp(files(k).name, ...
        '^sample(\d+)_(\w+?)_T([\d.]+)_tw([\d.]+)_p([\d.]+)\.mat$', 'tokens');
    if isempty(tok)
        fprintf('Skipping unrecognised file: %s\n', files(k).name);
        continue
    end
    sid = str2double(tok{1}{1});
    if ~isKey(t0_of, sid)
        fprintf('Skipping %s: no t0 known for sample %d\n', files(k).name, sid);
        continue
    end
    runs(end+1) = struct('file', files(k).name, 'sid', sid, 'geno', tok{1}{2}, ...
        'T', str2double(tok{1}{3}), 'tw', str2double(tok{1}{4}), ...
        'p', str2double(tok{1}{5}), 't0', t0_of(sid)); %#ok<SAGROW>
end

% per-run bounding boxes (mesh over t0..end, plus the simulated rod if asked)
box = nan(numel(runs), 6);       % xmin xmax ymin ymax zmin zmax
for r = 1:numel(runs)
    R  = runs(r);
    mb = load(sprintf('../data/meshes_spRS/meshes_sample%d.mat', R.sid));
    mb = mb.meshes;
    dref = load(sprintf('../data/refIV6_sample%d_walls3_cball.mat', R.sid), 'x');
    t_end = min(size(dref.x,1), numel(mb));
    b = [inf -inf inf -inf inf -inf];
    for t = R.t0:t_end
        v = mb{t}.v;
        b = [min(b(1),min(v(:,1))) max(b(2),max(v(:,1))) ...
             min(b(3),min(v(:,2))) max(b(4),max(v(:,2))) ...
             min(b(5),min(v(:,3))) max(b(6),max(v(:,3)))];
    end
    if fit_simulation
        S = load(fullfile(selected_dir, R.file), 'xsaves');
        L = sum(sqrt(sum(diff(squeeze(dref.x(R.t0,:,:)),1,1).^2, 2)));
        xs = S.xsaves * L;
        b = [min(b(1),min(min(xs(:,:,1)))) max(b(2),max(max(xs(:,:,1)))) ...
             min(b(3),min(min(xs(:,:,2)))) max(b(4),max(max(xs(:,:,2)))) ...
             min(b(5),min(min(xs(:,:,3)))) max(b(6),max(max(xs(:,:,3))))];
    end
    box(r,:) = b + axes_pad*[-1 1 -1 1 -1 1];
    clear mb
end
if shared_across_samples
    gbox = [min(box(:,1)) max(box(:,2)) min(box(:,3)) max(box(:,4)) ...
            min(box(:,5)) max(box(:,6))];
    box = repmat(gbox, numel(runs), 1);
    fprintf('shared axes limits: x [%.0f %.0f] y [%.0f %.0f] z [%.0f %.0f]\n', gbox);
end

views = {'ventral', 'side'};

% ---- render ----------------------------------------------------------------
fig = figure('Units','pixels', 'Position',[10 10 fig_px], 'Color','w', ...
             'PaperPositionMode','auto', 'Visible','off');

for r = 1:numel(runs)
    R = runs(r);
    if strcmpi(R.geno, 'WT')
        vik_c = vik(22,:);
    else
        vik_c = vik(234,:);
    end

    meshes = load(sprintf('../data/meshes_spRS/meshes_sample%d.mat', R.sid));
    meshes = meshes.meshes;
    data   = load(sprintf('../data/refIV6_sample%d_walls3_cball.mat', R.sid));
    data_x = data.x;  data_Q = data.Q;
    t_end  = min(size(data_x,1), numel(meshes));
    t0     = R.t0;

    S = load(fullfile(selected_dir, R.file));
    xsaves = S.xsaves;  Qsaves = S.Qsaves;

    % the n_saves snapshots span [0,T] <-> reference frames t0..num_t
    t_show = t0:t_end;
    sim_t  = linspace(t0, size(data_x,1), size(xsaves,1));
    L_t0   = sum(sqrt(sum(diff(squeeze(data_x(t0,:,:)),1,1).^2, 2)));
    xi = interp1(sim_t, squeeze(xsaves(:,:,1)), t_show) * L_t0;
    yi = interp1(sim_t, squeeze(xsaves(:,:,2)), t_show) * L_t0;
    zi = interp1(sim_t, squeeze(xsaves(:,:,3)), t_show) * L_t0;
    d1xi = interp1(sim_t, squeeze(Qsaves(:,:,director_idx,1)), t_show);
    d1yi = interp1(sim_t, squeeze(Qsaves(:,:,director_idx,2)), t_show);
    d1zi = interp1(sim_t, squeeze(Qsaves(:,:,director_idx,3)), t_show);

    out_dir = fullfile(out_root, erase(R.file, '.mat'));
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end
    fprintf('==== %s  (sample %d, %s, t0=%d, %d frames) -> %s\n', ...
            R.file, R.sid, R.geno, t0, numel(t_show), out_dir);

    for it = 1:numel(t_show)
        t = t_show(it);
        clf(fig);
        tiledlayout(fig, 1, numel(views), 'Padding','tight', 'TileSpacing','compact');

        for iv = 1:numel(views)
            nexttile;
            trisurf(triangulation(meshes{t}.f, meshes{t}.v), 'EdgeColor','none', ...
                    'FaceColor',0.9*[1 1 1], 'FaceAlpha',0.2);
            hold on

            [xx, yy, zz] = tubeplot(squeeze(data_x(t,:,:))', tube_r, tube_n);
            surf(xx, yy, zz, 'FaceColor',0.6*[1 1 1], 'EdgeColor','none');

            [xx, yy, zz] = tubeplot([xi(it,:); yi(it,:); zi(it,:)], tube_r, tube_n);
            surf(xx, yy, zz, 'FaceColor',vik_c, 'EdgeColor','none');

            if draw_quivers
                % every quiver3 argument forced to a column: data_x(t,sub,1)
                % is a 1xN row, d1(sub,1) an Nx1 column, and mixing them
                % errors with "Z and U must be the same size"
                if draw_ref_frame
                    d1r  = squeeze(data_Q(t, :, director_idx, :));
                    rsub = 1:quiver_stride:size(d1r, 1);
                    bx = data_x(t,rsub,1); by = data_x(t,rsub,2); bz = data_x(t,rsub,3);
                    quiver3(bx(:), by(:), bz(:), ...
                            d1r(rsub,1)*quiver_len, d1r(rsub,2)*quiver_len, ...
                            d1r(rsub,3)*quiver_len, 0, 'Color','k', 'LineWidth',quiver_lw);
                end
                ssub = 1:quiver_stride:size(d1xi, 2);
                sx = xi(it,ssub);   sy = yi(it,ssub);   sz = zi(it,ssub);
                ux = d1xi(it,ssub); uy = d1yi(it,ssub); uz = d1zi(it,ssub);
                quiver3(sx(:), sy(:), sz(:), ...
                        ux(:)*quiver_len, uy(:)*quiver_len, uz(:)*quiver_len, ...
                        0, 'Color',vik_c, 'LineWidth',quiver_lw);
            end

            axis(box(r,:));
            daspect([1 1 1]);
            axis off
            title(views{iv}, 'FontSize',9);
            lighting gouraud
            material dull
            if strcmp(views{iv}, 'side')
                view([0, 0]);
            else
                view(0, -90);
            end
            camlight('headlight');
        end

        % backslashes doubled: sprintf would turn \t into a tab and stop at the
        % invalid escape \e; the TeX interpreter then renders \eta\tau_E
        sgtitle(sprintf(['sample %d (%s)   T=%g, B_t/B=%g, \\eta\\tau_E=%g   ' ...
                         't = %g min   (frame %d)\ngrey = imaged, colour = simulation'], ...
                R.sid, R.geno, R.T, R.tw, R.p, (t - t0)*dt_min, t), 'FontSize',10);

        print(fig, fullfile(out_dir, sprintf('frame_%03d.png', it)), '-dpng', save_dpi);
    end
    clear meshes data S
end

close(fig)
fprintf('~done~\n');
