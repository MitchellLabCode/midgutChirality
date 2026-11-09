% byn[5] LVM Scoring Results
% Scoring of longitudinal visceral muscle (LVM) presence in byn^5 mutants.
close all
clear all

% ---- Data (read from Excel) ---------------------------------------------
% Expected layout (Sheet 1): column 1 = category name, column 2 = count.
% Rows are mutually exclusive categories, ordered by increasing LVM coverage
% (Absent -> Scattered cells in posterior -> Fibers extending to anterior).
inputFile = fullfile(fileparts(mfilename('fullpath')), 'utils', 'SUPPLEMENTAL_byn5_LVM_coverage.xlsx');
data = readtable(inputFile, 'Sheet', 1, 'VariableNamingRule', 'preserve');

rawLabels = string(data{:, 1})';
counts    = data{:, 2}';
assert(numel(counts) == 3, 'Expected 3 categories in %s, found %d.', inputFile, numel(counts));

% Word-wrap each label into a cell array of lines so every line is centered
% individually (a single string with \newline centers the block but
% left-aligns its lines).
maxLineChars = 12;
catLabels = cell(1, numel(rawLabels));
for k = 1:numel(rawLabels)
    words = split(rawLabels(k))';
    lines = {};
    cur = '';
    for w = words
        if isempty(cur)
            cur = char(w);
        elseif strlength(cur) + 1 + strlength(w) <= maxLineChars
            cur = [cur ' ' char(w)];
        else
            lines{end+1} = cur; %#ok<SAGROW>
            cur = char(w);
        end
    end
    lines{end+1} = cur; %#ok<SAGROW>
    catLabels{k} = lines;
end

nTotal = sum(counts);
assert(nTotal == 10, 'Category counts sum to %d, expected 10.', nTotal);

% Bar fill styles: 'solid' = black, 'open' = white with black outline,
% 'hatch' = white with black outline and black diagonal hatching.
fillStyle = {'solid', 'open', 'hatch'};

barWidth     = 0.6;    % data units
edgeWidth    = 0.75;   % pt, bar outline
hatchWidth   = 0.5;    % pt, hatch lines
hatchSpacing = 0.15;   % cm, perpendicular distance between hatch lines

% ---- Figure -------------------------------------------------------------
fig = figure('Color', 'w');
% R2025a+ applies the desktop (possibly dark) theme to figures; force light.
if ~isMATLABReleaseOlderThan("R2025a")
    theme(fig, "light");
end
set(fig, 'Units', 'centimeters', 'Position', [2, 2, 9, 9]);
set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [9, 9], 'PaperPosition', [0, 0, 9, 9]);
hold on;

yMax = 6;   % one unit of headroom above the tallest bar (5); ends on a tick
assert(yMax > max(counts), 'yMax (%d) must exceed the tallest bar (%d).', yMax, max(counts));

for k = 1:numel(counts)
    if strcmp(fillStyle{k}, 'solid')
        faceColor = 'k';
    else
        faceColor = 'w';
    end
    bar(k, counts(k), barWidth, 'FaceColor', faceColor, 'EdgeColor', 'k', ...
        'LineWidth', edgeWidth);
end

% ---- Formatting ---------------------------------------------------------
ax = gca;
set(ax, 'Box', 'off', 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
    'FontName', 'Arial', 'FontSize', 9, 'FontSizeMode', 'manual', ...
    'Layer', 'top', 'TickLength', [0.02, 0.02], ...
    'TickDir', 'in', 'TickDirMode', 'manual');   % ticks point into the plot only
set(ax, 'YTick', 0:1:yMax);
ylim([0, yMax]);
xlim([0.4, numel(counts) + 0.6]);

% No x tick marks: remove x ticks entirely and place category labels as
% text just below the x-axis.
set(ax, 'XTick', []);
labelOffset = 0.03 * yMax;   % gap between axis line and label (data units)
for k = 1:numel(catLabels)
    text(ax, k, -labelOffset, catLabels{k}, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'Interpreter', 'tex', 'FontName', 'Arial', 'FontSize', 9, ...
        'Color', 'k', 'Clipping', 'off');
end

% Fixed axes position so the text labels (which MATLAB doesn't reserve
% space for) aren't clipped at the figure's bottom edge. Units: cm.
% If the bottom label line is clipped on export, raise the 2nd value and
% lower the 4th by the same amount, e.g. [1.6, 2.5, 7.0, 5.5]. Bottom margin of 2.1 cm fits up to 4 label lines.
set(ax, 'Units', 'centimeters', 'Position', [1.6, 2.1, 7.0, 5.9]);
ylabel('Number of Embryos');

% ---- Hatching -----------------------------------------------------------
% MATLAB bars have no built-in hatch fill, so draw 45-degree lines clipped
% to each hatched bar. Lines are laid out in physical (cm) space so they are
% truly 45 degrees despite unequal x/y data scaling. This must run AFTER the
% axes position and limits are final; re-run it if you change either.
axPos = ax.Position;                  % [left bottom width height], cm
sx = axPos(3) / diff(ax.XLim);        % cm per x data unit
sy = axPos(4) / diff(ax.YLim);        % cm per y data unit
dc = hatchSpacing * sqrt(2);          % intercept step for 45-degree lines
for k = find(strcmp(fillStyle, 'hatch'))
    x0 = k - barWidth/2;
    W  = barWidth  * sx;              % bar width, cm
    H  = counts(k) * sy;              % bar height, cm
    % Lines v = u + c in bar-local cm coords (u right, v up), c in (-W, H)
    for c = (-W + dc/2):dc:H
        u1 = max(0, -c);
        u2 = min(W, H - c);
        if u2 > u1
            line(ax, x0 + [u1, u2]/sx, ([u1, u2] + c)/sy, ...
                'Color', 'k', 'LineWidth', hatchWidth);
        end
    end
    % Redraw outline on top so hatch line ends don't nibble the border
    rectangle(ax, 'Position', [x0, 0, barWidth, counts(k)], ...
        'EdgeColor', 'k', 'LineWidth', edgeWidth);
end

% Title: byn^5 italic (TeX \it ... \rm), not bold, Helvetica 9 pt.
ttl = title('\itbyn^{5}\rm LVM Coverage', ...
    'Interpreter', 'tex', 'FontName', 'Helvetica', 'FontWeight', 'normal');
ax.TitleFontSizeMultiplier = 1;   % default is 1.1x axes font size
ttl.FontSize = 9;

set(fig, 'Color', 'w');
hold off;

saveas(gcf, './figure_byn5LVMScoring.fig')
exportgraphics(gcf, './figure_byn5LVMScoring.png', 'Resolution', 400, 'BackgroundColor', 'white')
exportgraphics(gcf, './figure_byn5LVMScoring.pdf', 'ContentType', 'vector', 'BackgroundColor', 'white')