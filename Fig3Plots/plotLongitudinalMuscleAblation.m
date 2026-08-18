% Longitudinal Muscle Ablation - Midgut Inversion Chirality
close all
clear all
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'))

% Read data from the Excel file
inputFile = fullfile(fileparts(mfilename('fullpath')), 'utils', 'LongitudinalMuscleAblation.xlsx');
data = readtable(inputFile, 'Sheet', 1, 'VariableNamingRule', 'preserve');

% Extract conditions and data rows
conditions = strrep(data.Properties.VariableNames(2:end), '_', ' ');
wts   = data{strcmp(data.("Midgut Coiling"), 'Percent WT'), 2:end};
revs  = data{strcmp(data.("Midgut Coiling"), 'Percent inverted'), 2:end};
abns  = 1 - wts - revs;
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

% Fisher's exact test vs w[1118] (first condition)
wt_rev_count     = round(revs(1) * n_values(1));
wt_non_rev_count = n_values(1) - wt_rev_count;
pvals    = ones(numBars, 1);
sig_strs = repmat({''}, numBars, 1);

for i = 2:numBars
    cond_rev_count     = round(revs(i) * n_values(i));
    cond_non_rev_count = n_values(i) - cond_rev_count;
    tbl = [cond_rev_count, cond_non_rev_count;
           wt_rev_count,   wt_non_rev_count];
    [~, pval] = fishertest(tbl, 'Tail', 'both');
    pvals(i) = pval;
    if     pval <= 0.0001, sig_strs{i} = '****';
    elseif pval <= 0.001,  sig_strs{i} = '***';
    elseif pval <= 0.01,   sig_strs{i} = '**';
    elseif pval <= 0.05,   sig_strs{i} = '*';
    end
end

% Labels above bars
for i = 1:numBars
    rev_count = round(revs(i) * n_values(i));
    if isempty(sig_strs{i})
        label_str = sprintf('\n%d/%d', rev_count, n_values(i));
    else
        label_str = sprintf('%s\n%d/%d', sig_strs{i}, rev_count, n_values(i));
    end
    text(i, revs(i) + err_upper(i) + 0.05, label_str, ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
end

% Print p-value table
fprintf('\n%s\n', repmat('=', 80, 1));
fprintf('P-values for inverted proportions vs w[1118]:\n');
fprintf('%s\n', repmat('=', 80, 1));
fprintf('%-40s %15s %10s\n', 'Condition', 'p-value', 'Significance');
fprintf('%s\n', repmat('-', 80, 1));
fprintf('%-40s %15s %10s\n', conditions{1}, '(reference)', '');
for i = 2:numBars
    if pvals(i) < 0.001
        pval_str = sprintf('%.1e', pvals(i));
    else
        pval_str = sprintf('%.4f', pvals(i));
    end
    fprintf('%-40s %15s %10s\n', conditions{i}, pval_str, sig_strs{i});
end
fprintf('%s\n', repmat('=', 80, 1));
fprintf('Significance levels: * p<=0.05, ** p<=0.01, *** p<=0.001, **** p<=0.0001\n');

% Formatting
xticks(1:numBars);
xticklabels(conditions);
xtickangle(45);
title('Longitudinal Muscle Ablation');
ylabel('Fraction inverted');
ylim([0, 1.05]);
set(gca, 'YTick', [0, 0.25, 0.5, 0.75, 1.0]);
set(gca, 'Box', 'off', 'TickDir', 'in', 'XColor', 'k');
set(gcf, 'color', 'w');
hold off;
legend('inverted', 'abnormal', 'normal', 'location', 'northwest');
legend boxoff

saveas(gcf, './figure_longitudinalMuscleAblation.fig')
exportgraphics(gcf, './figure_longitudinalMuscleAblation.png', 'Resolution', 400)
