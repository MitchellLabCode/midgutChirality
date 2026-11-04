% byn[5] Live Scoring Results
close all
clear all
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'))

% Read data from the Excel file
inputFile = fullfile(fileparts(mfilename('fullpath')), 'utils', 'SUPPLEMENTAL_byn5livescoringresults.xlsx');
data = readtable(inputFile, 'Sheet', 1, 'VariableNamingRule', 'preserve');

devSuccess = data.("DevSuccess");
chirality  = data.("Chirality");
counts     = data.("NumEmbryos");

failCount     = sum(counts(strcmp(devSuccess, 'Fail to Form 3 Constrictions')));
normalCount   = sum(counts(strcmp(chirality, 'Normal')));
invertedCount = sum(counts(strcmp(chirality, 'Inverted')));
abnormalCount = sum(counts(strcmp(chirality, 'Abnormal')));
formedCount   = normalCount + invertedCount + abnormalCount;

% Build bar data: two bars at x=1 and x=2
% Bar 1 (Fail): solid gray, one segment
% Bar 2 (Formed): stacked as abnormal (bottom), inverted (middle), normal (top)
barData = [failCount, 0,             0,             0; ...
           0,         abnormalCount,  invertedCount, normalCount];

fig = figure;
set(fig, 'Units', 'centimeters', 'Position', [2, 2, 9, 9]);
set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [9, 9], 'PaperPosition', [0, 0, 9, 9]);
hold on;

% Draw "Fail" bar manually as a single-color bar
failBar = bar(1, failCount, 0.6, 'FaceColor', [55/255, 55/255, 55/255], ...
    'EdgeColor', 'none');

% Draw "Formed" bar as stacked
formedData = [invertedCount, abnormalCount, normalCount];
bFormed = bar(2, formedData, 0.6, 'stacked', 'EdgeColor', 'none');
bFormed(1).FaceColor = [188/255, 64/255,  58/255];  % inverted (bottom)
bFormed(2).FaceColor = [243/255, 172/255, 172/255]; % abnormal (middle)
bFormed(3).FaceColor = [77/255,  106/255, 164/255]; % normal (top)

% Formatting
set(gca, 'Box', 'off', 'TickDir', 'in', 'XColor', 'k', 'YColor', 'k', ...
    'FontName', 'Arial', 'FontSize', 9, 'FontSizeMode', 'manual');
set(gca, 'XTick', [1, 2]);
set(gca, 'XTickLabel', {'Fail to Form 3 Constrictions', 'Scorable'});
ylim([0, max(failCount, formedCount) + 3]);
xlim([0.4, 2.6]);
ylabel('Number of Embryos');
title('byn^{5} Live Scoring Results');
set(gcf, 'color', 'w');
lg = legend([bFormed(3), bFormed(2), bFormed(1)], {'Normal', 'Abnormal', 'Inverted'}, ...
    'Location', 'eastoutside');
lg.FontSize = 9;
lg.FontName = 'Arial';
legend boxoff
hold off;

saveas(gcf, './figure_byn5LiveScoring.fig')
exportgraphics(gcf, './figure_byn5LiveScoring.png', 'Resolution', 400)
exportgraphics(gcf, './figure_byn5LiveScoring.pdf', 'ContentType', 'vector')
