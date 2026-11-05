function ax = colorbar_extend(fig, pos, cmap, cl, ticks, label)
% colorbar_extend  Vertical colorbar with triangular end caps, like
% matplotlib's colorbar(extend='both'): the caps mark that values beyond
% the limits are shown in the end colors (saturated).
%
%   ax = colorbar_extend(fig, pos, cmap, cl, ticks, label)
%
% fig   : figure handle
% pos   : [left bottom width height] of the bar body, normalized figure units
% cmap  : n x 3 colormap
% cl    : [cmin cmax] color limits
% ticks : tick values (within cl)
% label : colorbar label (TeX), '' for none

h = 0.12 * diff(cl);                   % cap height, in data units
ax = axes(fig, 'Position', [pos(1), pos(2) - 0.12*pos(4), pos(3), 1.24*pos(4)]);
hold(ax, 'on')

% body: the colormap as an RGB image, so it does not depend on the axes colormap
% (XData/YData are pixel centres: one column centred at 0.5 spans 0..1, and
% n rows centred half a step in from the limits span exactly cl)
n = size(cmap, 1);
img = reshape(cmap, [n, 1, 3]);
dy = diff(cl) / n;
image(ax, 'XData', 0.5, 'YData', [cl(1) + dy/2, cl(2) - dy/2], 'CData', img);

% end caps in the end colors
patch(ax, [0 1 0.5], [cl(2) cl(2) cl(2)+h], cmap(end,:), 'EdgeColor', 'k', 'LineWidth', 0.5);
patch(ax, [0 1 0.5], [cl(1) cl(1) cl(1)-h], cmap(1,:),   'EdgeColor', 'k', 'LineWidth', 0.5);
% outline of the body
plot(ax, [0 1 1 0 0], [cl(1) cl(1) cl(2) cl(2) cl(1)], 'k', 'LineWidth', 0.5);

% ticks and labels on the right
for k = 1:numel(ticks)
    plot(ax, [1 1.25], ticks(k)*[1 1], 'k', 'LineWidth', 0.5);
    text(ax, 1.45, ticks(k), num2str(ticks(k)), 'VerticalAlignment', 'middle', ...
         'FontSize', 10);
end
if ~isempty(label)
    text(ax, 3.2, mean(cl), label, 'Rotation', 90, 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', 11);
end

xlim(ax, [0 1]); ylim(ax, [cl(1)-h, cl(2)+h]);
set(ax, 'Visible', 'off', 'Clipping', 'off', 'YDir', 'normal')
end
