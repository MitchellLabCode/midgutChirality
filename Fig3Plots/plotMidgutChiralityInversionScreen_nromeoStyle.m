% Midgut Inversion Chirality Myo1C overexpression Myo1D KD_KO
% Plot results from MidgutInversionChirality_Myo1Coverexpression_Myo1DKDKO.xlsx
close all
clear all
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'))

% Read data from the Excel file
inputFile = fullfile(fileparts(mfilename('fullpath')), 'utils', 'MidgutInversionChirality_Myo1Coverexpression_Myo1DKDKO.xlsx');
sheet = 1; % Adjust if needed
data = readtable(inputFile, 'Sheet', sheet, 'VariableNamingRule', 'preserve');


% Extract relevant data
expt = data.Properties.VariableNames(2:end); % Extract condition names (skip first column)
wt = data{strcmp(data.("Midgut Coiling"), 'Percent WT'), 2:end}; % WT proportion
rev = data{strcmp(data.("Midgut Coiling"), 'Percent inverted'), 2:end}; % Inverted proportion
abnormal = 1 - wt - rev ;
n_values = data{strcmp(data.("Midgut Coiling"), 'n'), 2:end}; % Sample sizes

% Number of conditions
dmyk = 1 ;
conditions = {} ;
wts = [] ;
revs = [] ;
anormal = [] ;
n_vals = [] ;
for cID = 1:length(expt)
    if isempty(strfind(expt{cID}, 'ontrol'))
        conditions{dmyk} = strrep(expt{cID}, '_', ' ') ;
        wts(dmyk) = wt(cID) ;
        revs(dmyk) = rev(cID) ;
        abns(dmyk) = abnormal(cID) ;
        n_vals(dmyk) = n_values(cID) ;
        dmyk = dmyk + 1 ;
    end
end
numBars = length(conditions);

% Reorder the conditions
order = [1,3,5:9, 11:14,16:17, 2, 4, 19,20] ;
conditions = conditions(order);
wts = wts(order);
revs = revs(order);
abns = abns(order);
n_values = n_vals(order);
numBars = length(conditions); % Recalculate after reordering


% Create figure
figure; hold on;
b = bar([revs ; abns]', 'stacked', 'EdgeColor','none');

% Set colors (red for inverted on bottom, gray for WT on top)
b(1).FaceColor = [204/255, 51/255, 51/255];%[0., 0.52, 0.52]; %[0.75, 0, 0]; % Inverted (bottom)
b(2).FaceColor = [0.8, 0.8, 0.8]; % Abnormal (middle)
%b(3).FaceColor = [0.5, 0.5, 0.5]; % WT (top)

% Calculate and add error bars for red (inverted) portions
err_lower = zeros(numBars, 1);
err_upper = zeros(numBars, 1);

for i = 1:numBars
    % Convert proportions to counts
    rev_count = round(revs(i) * n_values(i));
    abn_count = round(abns(i) * n_values(i));
    wt_count = round(wts(i) * n_values(i));

    % Get uncertainty estimates
    if (rev_count + abn_count + wt_count) > 0
        stats = threeway_uncertainty_wilson(rev_count, abn_count, wt_count);
        % Distance DOWN from revs to lower CI, and UP from revs to upper CI
        err_lower(i) = revs(i) - stats.ci95(1, 1);  % Distance to lower bound
        err_upper(i) = stats.ci95(1, 2) - revs(i);  % Distance to upper bound
    end
end

% Add error bars to the red (inverted) portion - centered on the red bar height
x_pos = 1:numBars;
errorbar(x_pos, revs, err_lower, err_upper, 'k', 'LineStyle', 'none', ...
    'CapSize', 5, 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Calculate p-values for inverted (red) proportions vs WT (first condition)
pvals = zeros(numBars, 1);
sig_strs = cell(numBars, 1);

% Get WT (first condition) counts
wt_rev_count = round(revs(1) * n_values(1));
wt_non_rev_count = n_values(1) - wt_rev_count;

% WT is reference, no p-value
pvals(1) = 1;
sig_strs{1} = '';

% Compare each condition to WT
for i = 2:numBars
    % Convert proportions to counts for condition i
    cond_rev_count = round(revs(i) * n_values(i));
    cond_non_rev_count = n_values(i) - cond_rev_count;

    % 2x2 contingency table: [condition_i_inverted, condition_i_non_inverted;
    %                         wt_inverted, wt_non_inverted]
    table = [cond_rev_count, cond_non_rev_count;
             wt_rev_count, wt_non_rev_count];

    % Fisher's exact test (two-tailed)
    [~, pval] = fishertest(table, 'Tail', 'both');
    pvals(i) = pval;

    % Assign significance asterisks
    if pval <= 0.0001
        sig_strs{i} = '****';
    elseif pval <= 0.001
        sig_strs{i} = '***';
    elseif pval <= 0.01
        sig_strs{i} = '**';
    elseif pval <= 0.05
        sig_strs{i} = '*';
    else
        sig_strs{i} = '';
    end
end

% Display sample size and significance above each bar
for i = 1:numBars
    rev_count = round(revs(i) * n_values(i));
    if isempty(sig_strs{i})
        label_str = sprintf('\n%d/%d',rev_count, n_values(i));
    else
        label_str = sprintf('%s\n%d/%d', sig_strs{i}, rev_count, n_values(i));
    end
    text(i, revs(i)+err_upper(i)+0.05, label_str, 'HorizontalAlignment', 'center', 'FontSize', 9);
end


% Report p-values on command line
fprintf('\n%s\n', repmat('=', 80, 1));
fprintf('P-values for inverted proportions vs WT (first condition):\n');
fprintf('%s\n', repmat('=', 80, 1));
fprintf('%-40s %15s %10s\n', 'Condition', 'p-value', 'Significance');
fprintf('%s\n', repmat('-', 80, 1));

fprintf('%-40s %15s %10s\n', conditions{1}, '(reference)', '');

for i = 2:numBars
    pval = pvals(i);
    % Format p-value: use scientific notation for very small values
    if pval < 0.001
        pval_str = sprintf('%.1e', pval);
    else
        pval_str = sprintf('%.4f', pval);
    end

    fprintf('%-40s %15s %10s\n', conditions{i}, pval_str, sig_strs{i});
end
fprintf('%s\n', repmat('=', 80, 1));
fprintf('Significance levels: * p<=0.05, ** p<=0.01, *** p<=0.001, **** p<=0.0001, ns = not significant\n');

% Formatting
xticks(1:numBars);
xticklabels(conditions);
xtickangle(45); % Rotate x-axis labels for readability
% xlabel('Condition');

% ylabel('Proportion');
title('Midgut Inversion');
ylim([0, 1.2]); % Adjust to fit labels

set(gca, 'YTick', [0, 0.25, 0.5, 0.75, 1.0]); % Removes y-axis ticks
ylim([0., 1.05])
ylabel('Fraction inverted')
%set(gca, 'YColor', 'none'); % Turns off the y-axis line
set(gca, 'Box', 'off', 'TickDir', 'in', 'XColor', 'k')%, 'YColor', 'none')
set(gcf, 'color', 'w')
hold off;
legend('inverted', 'abnormal', 'normal', 'location', 'northwest') %'eastoutside')
legend boxoff
saveas(gcf, './figure_midgutChiralityInversionMyo1CMyo1D.fig')
exportgraphics(gcf,'./figure_midgutChiralityInversionMyo1CMyo1D.png','Resolution',400)
