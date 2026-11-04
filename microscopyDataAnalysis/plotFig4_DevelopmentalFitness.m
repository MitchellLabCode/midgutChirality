% Fig 4: Developmental Fitness — Hatching, Pupation, Eclosion
% Wilson 95% CIs, two-sided Fisher's exact tests
close all
clear all
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'))

inputFile = fullfile(fileparts(mfilename('fullpath')), 'utils', 'DevelopmentalFitness_Fig4.xlsx');

% --- Read data ---
hatchData  = readcell(inputFile, 'Range', 'A3:C4');
pupEclData = readcell(inputFile, 'Range', 'A9:D10');

groups = {'Normal'; 'Inverted'};
stages = {'Hatching', 'Pupation', 'Eclosion'};

% rows = [Normal; Inverted], cols = [Hatching, Pupation, Eclosion]
successes = [hatchData{1,3}, pupEclData{1,3}, pupEclData{1,4};
             hatchData{2,3}, pupEclData{2,3}, pupEclData{2,4}];
trials    = [hatchData{1,2}, pupEclData{1,2}, pupEclData{1,3};
             hatchData{2,2}, pupEclData{2,2}, pupEclData{2,3}];

alpha = 0.05;
nGroups = 2;
nStages = 3;

% --- Compute CIs ---
pHat = zeros(nGroups, nStages);
wLo  = zeros(nGroups, nStages);  wHi  = zeros(nGroups, nStages);

for g = 1:nGroups
    for s = 1:nStages
        [pHat(g,s), wLo(g,s), wHi(g,s)] = wilsonCI(successes(g,s), trials(g,s), alpha);
    end
end

% --- Print CI table ---
Population = {}; Stage = {}; Succ = []; Tri = [];
pH = []; wL = []; wH = [];
for g = 1:nGroups
    for s = 1:nStages
        Population{end+1,1} = groups{g};
        Stage{end+1,1}      = stages{s};
        Succ(end+1,1) = successes(g,s);
        Tri(end+1,1)  = trials(g,s);
        pH(end+1,1)   = pHat(g,s);
        wL(end+1,1)   = wLo(g,s);
        wH(end+1,1)   = wHi(g,s);
    end
end
T = table(string(Population), string(Stage), Succ, Tri, pH, wL, wH, ...
    'VariableNames', {'Population','Stage','Successes','Trials', ...
    'pHat','Wilson_lo','Wilson_hi'});
disp(T);

% --- Two-sided Fisher's exact tests (Normal vs Inverted at each stage) ---
fprintf('\n%s\n', repmat('=', 70, 1));
fprintf('Two-sided Fisher''s exact test: Normal vs Inverted\n');
fprintf('%s\n', repmat('=', 70, 1));
fprintf('%-12s %10s %10s\n', 'Stage', 'p-value', 'Sig');
fprintf('%s\n', repmat('-', 70, 1));
pvals = zeros(1, nStages);
sigs  = cell(1, nStages);
for s = 1:nStages
    tbl = [successes(1,s), trials(1,s) - successes(1,s);
           successes(2,s), trials(2,s) - successes(2,s)];
    [~, pvals(s)] = fishertest(tbl, 'Tail', 'both');
    if     pvals(s) <= 0.0001, sigs{s} = '****';
    elseif pvals(s) <= 0.001,  sigs{s} = '***';
    elseif pvals(s) <= 0.01,   sigs{s} = '**';
    elseif pvals(s) <= 0.05,   sigs{s} = '*';
    else,                      sigs{s} = 'ns';
    end
    if pvals(s) < 0.001, pv_str = sprintf('%.2e', pvals(s));
    else,                pv_str = sprintf('%.4f', pvals(s));
    end
    fprintf('%-12s %10s %10s\n', stages{s}, pv_str, sigs{s});
end
fprintf('%s\n', repmat('=', 70, 1));
fprintf('Significance: ns p>0.05, * p<=0.05, ** p<=0.01, *** p<=0.001, **** p<=0.0001\n');

% --- Butterfly bar plot with Wilson 95% CIs ---
% Normal extends left from center, Inverted extends right
norm_p  = pHat(1,:) * 100;
norm_lo = wLo(1,:)  * 100;
norm_hi = wHi(1,:)  * 100;
inv_p   = pHat(2,:) * 100;
inv_lo  = wLo(2,:)  * 100;
inv_hi  = wHi(2,:)  * 100;

% Y positions: Hatching=3 (top), Pupation=2, Eclosion=1 (bottom)
yPos = [3, 2, 1];
barH = 0.5;

cNorm = [68, 107, 169] / 255;
cInv  = [188,  64,  58] / 255;

figure;
hold on;

for s = 1:nStages
    % Normal bars (negative direction = left)
    barh(yPos(s), -norm_p(s), barH, 'FaceColor', cNorm, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    % Inverted bars (positive direction = right)
    barh(yPos(s), inv_p(s),   barH, 'FaceColor', cInv,  'EdgeColor', 'none', 'HandleVisibility', 'off');

    % Error bars — Normal (left side)
    errorbar(-norm_p(s), yPos(s), [], [], norm_hi(s) - norm_p(s), norm_p(s) - norm_lo(s), ...
        'k', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 6, 'HandleVisibility', 'off');
    % Error bars — Inverted (right side)
    errorbar(inv_p(s), yPos(s), [], [], inv_p(s) - inv_lo(s), inv_hi(s) - inv_p(s), ...
        'k', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 6, 'HandleVisibility', 'off');

    % Significance label to the right
    text(max(inv_hi(s), inv_p(s)) + 3, yPos(s), sigs{s}, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 9);
end

% Invisible patches for legend
patch(NaN, NaN, cNorm, 'EdgeColor', 'none', 'DisplayName', 'Normal');
patch(NaN, NaN, cInv,  'EdgeColor', 'none', 'DisplayName', 'Inverted');

% Center line
plot([0 0], [0.3 3.7], 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');

% Column headers
text(-50, 3.5, 'Normal',   'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10);
text( 50, 3.5, 'Inverted', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10);

hold off;

% X-axis: 100 on left, 0 in center, 100 on right
set(gca, 'XTick', [-100, -50, 0, 50, 100], ...
         'XTickLabel', {'100', '50', '0', '50', '100'});
xlim([-110, 120]);
xlabel('Success rate (%)');

% Y-axis: stage labels
set(gca, 'YTick', [1, 2, 3], 'YTickLabel', {'Eclosion', 'Pupation', 'Hatching'});
ylim([0.3, 3.7]);

set(gca, 'Box', 'on', 'TickDir', 'in');
set(gcf, 'color', 'w');
lg = legend('Normal', 'Inverted', 'Location', 'southoutside', 'Orientation', 'horizontal');
legend boxoff

outDir = fileparts(mfilename('fullpath'));
saveas(gcf, fullfile(outDir, 'figure_DevelopmentalFitness_Fig4.fig'))
exportgraphics(gcf, fullfile(outDir, 'figure_DevelopmentalFitness_Fig4.png'), 'Resolution', 400)
exportgraphics(gcf, fullfile(outDir, 'figure_DevelopmentalFitness_Fig4.pdf'), 'ContentType', 'vector')
