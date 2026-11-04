% Midgut Chirality - byn[5] live scoring vs w[1118]
% Statistical test: Fisher's exact test (two-sided), single comparison
close all
clear all
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'))

% Read data from the Excel file
inputFile = fullfile(fileparts(mfilename('fullpath')), 'utils', 'MidgutChirality_byn5LiveScoringVsWT.xlsx');
data = readtable(inputFile, 'Sheet', 1, 'VariableNamingRule', 'preserve'); 

% Extract conditions and data rows
conditions = strrep(data.Properties.VariableNames(2:end), '_', ' ');
wts      = data{strcmp(data.("Midgut Coiling"), 'Percent WT'), 2:end};
revs     = data{strcmp(data.("Midgut Coiling"), 'Percent inverted'), 2:end};
abns     = 1 - wts - revs;
n_values = data{strcmp(data.("Midgut Coiling"), 'n'), 2:end};

numBars = length(conditions);

% Create figure
figure; hold on;
b = bar([revs; abns]', 'stacked', 'EdgeColor', 'none');
b(1).FaceColor = [204/255, 51/255, 51/255]; % inverted
b(2).FaceColor = [0.8, 0.8, 0.8];           % abnormal

% Wilson 95% CI error bars on inverted proportion
err_lower = zeros(numBars, 1);
err_upper = zeros(numBars, 1);
for i = 1:numBars
    rev_count = round(revs(i) * n_values(i));
    abn_count = round(abns(i) * n_values(i));
    wt_count  = round(wts(i)  * n_values(i));
    if (rev_count + abn_count + wt_count) > 0
        stats = threeway_uncertainty_wilson(rev_count, abn_count, wt_count);
        err_lower(i) = revs(i) - stats.ci95(1, 1);
        err_upper(i) = stats.ci95(1, 2) - revs(i);
    end
end
errorbar(1:numBars, revs, err_lower, err_upper, 'k', 'LineStyle', 'none', ...
    'CapSize', 5, 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Fisher's exact test: w[1118] vs byn[5] (two-sided, single comparison)
c1 = round(revs(1) * n_values(1));
c2 = round(revs(2) * n_values(2));
tbl = [c1, n_values(1)-c1; c2, n_values(2)-c2];
[~, pval] = fishertest(tbl, 'Tail', 'both');

if     pval <= 0.0001, sig = '****';
elseif pval <= 0.001,  sig = '***';
elseif pval <= 0.01,   sig = '**';
elseif pval <= 0.05,   sig = '*';
else,                  sig = 'ns';
end

% n/count labels above each bar
for i = 1:numBars
    rev_count = round(revs(i) * n_values(i));
    text(i, revs(i) + err_upper(i) + 0.02, sprintf('%d/%d', rev_count, n_values(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 8);
end

% Draw significance bracket
tickLen    = 0.015;
baseHeight = max(revs(:) + err_upper(:)) + 0.08;
y = baseHeight;
plot([1 2], [y y],         'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
plot([1 1], [y-tickLen y], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
plot([2 2], [y-tickLen y], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
text(1.5, y+0.01, sig, 'HorizontalAlignment', 'center', 'FontSize', 9);
needed_ylim = y + 0.08;

% Print p-value table
fprintf('\n%s\n', repmat('=', 60, 1));
fprintf('Fisher''s exact test (two-sided): %s vs %s\n', conditions{1}, conditions{2});
fprintf('%s\n', repmat('=', 60, 1));
if pval < 0.001, pv_str = sprintf('%.2e', pval); else, pv_str = sprintf('%.4f', pval); end
fprintf('p-value: %s  %s\n', pv_str, sig);
fprintf('%s\n', repmat('=', 60, 1));
fprintf('Significance: ns p>0.05, * p<=0.05, ** p<=0.01, *** p<=0.001, **** p<=0.0001\n');

% Formatting
xticks(1:numBars);
xticklabels(conditions);
xtickangle(45);
title('Midgut Chirality byn^5');
ylabel('Fraction inverted');
ylim([0, needed_ylim]);
set(gca, 'YTick', [0, 0.25, 0.5, 0.75, 1.0]);
set(gca, 'Box', 'off', 'TickDir', 'in', 'XColor', 'k');
set(gcf, 'color', 'w');
hold off;
legend('inverted', 'abnormal', 'location', 'northwest');
legend boxoff

saveas(gcf, './figure_midgutChirality_byn5.fig')
exportgraphics(gcf, './figure_midgutChirality_byn5.png', 'Resolution', 400)
exportgraphics(gcf, './figure_midgutChirality_byn5.pdf', 'ContentType', 'vector')
