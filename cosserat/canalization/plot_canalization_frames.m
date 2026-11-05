% Frame-by-frame PNGs of a canalization sweep: selected alpha runs evolving
% side by side, plus a writhe-vs-length panel that grows with time.
%
% Layout of every frame:
%   left  : one column per selected alpha, ventral (top) and side (bottom) tube
%           views at the same snapshot, coloured by alpha
%   right : Wr vs L/L0 for every run (minus the discontinuous one), each curve
%           drawn up to the current snapshot; the selected runs carry a marker
%           at their current point. Axes as in the paper figure: Wr along the
%           top, L/L0 increasing downward, a line at Wr = 0.
%
% Writhe comes from precompute_writhe.py (formula of ../src/measurements.py,
% verified against measurements.compute_wr), which also flags the run with a
% writhe discontinuity. Run that first:
%     python precompute_writhe.py 202609022_canalization_box2f_s6s1
%
% Colours: alpha = 1 - 2*tau on vik_r over [-1.2, 1.2].
% Positions are in units of the initial length (L0 = 1), as in the simulation.
%
% Output: <folder>/frames_wr/frame_NNN.png. Run from canalization/.

addpath(".")

% ======================= CONFIG ============================================
folder_name = '202609022_canalization_box2f_s6s1';
alpha_show  = [1, 0.5, 0, -0.5, -1];   % snapped to the nearest available run

views = {'ventral', 'side'};

tube_r = 0.013;                  % in L0 units (~5 um for L0 ~ 380 um)
tube_n = 24;
axes_pad = 0.03;                 % in L0 units

draw_quivers  = true;
director_idx  = 1;
quiver_len    = 0.065;           % in L0 units (~25 um); fixed, AutoScale off
quiver_stride = 6;
quiver_lw     = 1.0;

wr_lw     = 2.5;                 % writhe curve width
wr_xticks = [-0.5 0 0.5];
fig_px    = [1700 700];          % fixed pixel size -> every PNG identical size
save_dpi  = '-r150';

% ---- layout, normalized figure units ----
% The tube axes are sized to their content: daspect is 1:1 and the rods are
% ~2.2x longer than tall, so in an even tile grid ~70% of each tube tile was
% blank above and below the rod, which pushed the two view rows far apart.
% Sizing the axes to the content puts the rows exactly tube_vgap apart.
tube_region = [0.02 0.06 0.60 0.82];   % [left bottom width height] for the tubes
tube_hgap   = 0.006;                    % between alpha columns
tube_vgap   = 0.04;                     % between the ventral and side rows
wr_pos      = [0.69 0.10 0.28 0.72];   % writhe panel (inner position)
% ===========================================================================

W = load(fullfile(folder_name, 'wr_vs_len.mat'));
alpha    = W.alpha(:);
wr       = W.wr;                 % n_runs x n_frames
lr       = W.len_ratio;          % n_runs x n_frames, L/L0
snap_idx = W.snap_idx(:)';       % 1-based snapshot of each frame
excluded = logical(W.excluded(:));
files    = cellstr(W.files);
n_frames = numel(snap_idx);

vik = load('vik.mat').vik;
nv  = size(vik, 1);
% Normalize(-1.2, 1.2) on vik_r, i.e. vik read backwards
col_of = @(a) vik(1 + round((1 - (min(max(a,-1.2),1.2) + 1.2)/2.4) * (nv-1)), :);

% ---- selected runs: snap to the nearest run that is NOT excluded ----------
cand = find(~excluded);
sel  = zeros(size(alpha_show));
for k = 1:numel(alpha_show)
    [~, m] = min(abs(alpha(cand) - alpha_show(k)));
    sel(k) = cand(m);
    if abs(alpha(sel(k)) - alpha_show(k)) > 1e-9
        fprintf('alpha %g -> nearest available %.3f\n', alpha_show(k), alpha(sel(k)));
    end
end
sel = unique(sel, 'stable');
n_sel = numel(sel);
if any(excluded)
    fprintf('excluded from the writhe panel: alpha = %s (writhe discontinuity)\n', ...
            mat2str(round(alpha(excluded)', 3)));
end

% ---- load the selected runs; one shared tube box over all of them ---------
X = cell(n_sel,1);  Q = cell(n_sel,1);
b = [inf -inf inf -inf inf -inf];
for k = 1:n_sel
    S = load(fullfile(folder_name, files{sel(k)}));
    X{k} = S.xsaves;  Q{k} = S.Qsaves;
    xs = S.xsaves(snap_idx, :, :);
    b = [min(b(1),min(min(xs(:,:,1)))) max(b(2),max(max(xs(:,:,1)))) ...
         min(b(3),min(min(xs(:,:,2)))) max(b(4),max(max(xs(:,:,2)))) ...
         min(b(5),min(min(xs(:,:,3)))) max(b(6),max(max(xs(:,:,3))))];
end
tube_box = b + axes_pad*[-1 1 -1 1 -1 1];

% ---- tube axes positions, sized to the content aspect ---------------------
R      = tube_box([2 4 6]) - tube_box([1 3 5]);     % x, y, z ranges
asp    = fig_px(1) / fig_px(2);                     % normalized width -> height
view_h = zeros(1, numel(views));                    % content height/width per view
for iv = 1:numel(views)
    if strcmp(views{iv}, 'side')
        view_h(iv) = R(3) / R(1);                   % side: x across, z up
    else
        view_h(iv) = R(2) / R(1);                   % ventral/dorsal: x across, y up
    end
end
nv_rows = numel(views);
cw = (tube_region(3) - (n_sel-1)*tube_hgap) / n_sel;          % column width
blk = cw*asp*sum(view_h) + (nv_rows-1)*tube_vgap;            % block height
if blk > tube_region(4)                                       % shrink to fit
    cw  = (tube_region(4) - (nv_rows-1)*tube_vgap) / (asp*sum(view_h));
    blk = tube_region(4);
end
x0 = tube_region(1) + (tube_region(3) - (n_sel*cw + (n_sel-1)*tube_hgap)) / 2;
y  = tube_region(2) + (tube_region(4) + blk) / 2;            % top of the block
tube_pos = zeros(nv_rows, n_sel, 4);
for iv = 1:nv_rows
    h = cw*asp*view_h(iv);
    y = y - h;
    for k = 1:n_sel
        tube_pos(iv,k,:) = [x0 + (k-1)*(cw + tube_hgap), y, cw, h];
    end
    y = y - tube_vgap;
end

% ---- fixed writhe-panel limits, over every plotted (non-excluded) curve ---
keep = ~excluded;
wx = wr(keep,:);  wy = lr(keep,:);
wr_xlim = [min(wx(:)) max(wx(:))] + 0.05*range(wx(:))*[-1 1];
wr_ylim = [min(wy(:)) max(wy(:))] + 0.03*range(wy(:))*[-1 1];
[~, draw_order] = sort(alpha);              % draw curves in alpha order

out_dir = fullfile(folder_name, 'frames_wr');
if ~isfolder(out_dir)
    mkdir(out_dir);
end

fig = figure('Units','pixels', 'Position',[10 10 fig_px], 'Color','w', ...
             'PaperPositionMode','auto', 'Visible','off');

for f = 1:n_frames
    s = snap_idx(f);
    clf(fig);

    % ---- tubes: rows = views, cols = selected alpha -----------------------
    for iv = 1:numel(views)
        for k = 1:n_sel
            axes(fig, 'Position', squeeze(tube_pos(iv,k,:))');
            c  = col_of(alpha(sel(k)));
            xk = squeeze(X{k}(s,:,:))';             % 3 x N
            [xx, yy, zz] = tubeplot(xk, tube_r, tube_n);
            surf(xx, yy, zz, 'FaceColor',c, 'EdgeColor','none');
            hold on
            if draw_quivers
                ssub = 1:quiver_stride:size(Q{k}, 2);
                d1   = squeeze(Q{k}(s, ssub, director_idx, :));
                bx = X{k}(s,ssub,1); by = X{k}(s,ssub,2); bz = X{k}(s,ssub,3);
                quiver3(bx(:), by(:), bz(:), d1(:,1)*quiver_len, ...
                        d1(:,2)*quiver_len, d1(:,3)*quiver_len, 0, ...
                        'Color',0.6*c, 'LineWidth',quiver_lw);
            end
            axis(tube_box);
            daspect([1 1 1]);
            axis off
            lighting gouraud
            material dull
            if strcmp(views{iv}, 'side')
                view([0, 0]);
            elseif strcmp(views{iv}, 'dorsal')
                view(0, 90);
            else
                view(0, -90);
            end
            camlight('headlight');
            if iv == 1
                title(sprintf('\\alpha = %.2f', alpha(sel(k))), 'Color','k', 'FontSize',10);
            elseif k == 1
                title(views{iv}, 'FontSize',9);
            end
        end
    end

    % ---- writhe vs length, grown to the current frame ----------------------
    ax = axes(fig, 'Position', wr_pos);
    hold(ax, 'on');
    xline(ax, 0, 'Color',0.5*[1 1 1], 'LineWidth',1);
    for r = draw_order(:)'
        if excluded(r)
            continue
        end
        plot(ax, wr(r,1:f), lr(r,1:f), '-', 'Color',col_of(alpha(r)), 'LineWidth',wr_lw);
    end
    for k = 1:n_sel
        plot(ax, wr(sel(k),f), lr(sel(k),f), 'o', 'MarkerSize',8, ...
             'MarkerFaceColor',col_of(alpha(sel(k))), 'MarkerEdgeColor','k');
    end
    set(ax, 'YDir','reverse', 'XAxisLocation','top', 'TickDir','in', ...
            'Box','on', 'FontSize',12, 'XTick',wr_xticks, 'LineWidth',1);
    xlim(ax, wr_xlim);  ylim(ax, wr_ylim);
    xlabel(ax, 'Writhe Wr');
    ylabel(ax, 'Length L/L_0');
    axis(ax, 'square');

    sgtitle(fig, sprintf('%s   snapshot %d / %d   median L/L_0 = %.2f', ...
            strrep(folder_name,'_','\_'), s, size(X{1},1), median(lr(keep,f))), ...
            'FontSize',10);

    print(fig, fullfile(out_dir, sprintf('frame_%03d.png', f)), '-dpng', save_dpi);
end

close(fig)
fprintf('wrote %d frames to %s\n', n_frames, out_dir);
