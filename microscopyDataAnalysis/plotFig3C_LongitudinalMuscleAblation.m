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

% Specified pairwise Fisher's exact tests (two-sided) with Holm-Bonferroni correction
defined_pairs = [1,2; 1,3; 2,3; 3,4; 4,5; 1,5; 1,4];
num_pairs = size(defined_pairs, 1);

pair_list = defined_pairs;
pair_sig  = repmat({''}, num_pairs, 1);
raw_pval  = zeros(num_pairs, 1);
pair_pval = ones(num_pairs, 1);

for k = 1:num_pairs
    i = defined_pairs(k, 1);
    j = defined_pairs(k, 2);
    ci = round(revs(i) * n_values(i));
    cj = round(revs(j) * n_values(j));
    tbl = [ci, n_values(i)-ci; cj, n_values(j)-cj];
    [~, pval] = fishertest(tbl, 'Tail', 'both');
    raw_pval(k) = pval;
end

[sorted_p, sort_idx] = sort(raw_pval);
running_max = 0;
adj = zeros(num_pairs, 1);
for r = 1:num_pairs
    corrected   = min(1, sorted_p(r) * (num_pairs - r + 1));
    running_max = max(running_max, corrected);
    adj(r)      = running_max;
end
pair_pval(sort_idx) = adj;

for k = 1:num_pairs
    if     pair_pval(k) <= 0.0001, pair_sig{k} = '****';
    elseif pair_pval(k) <= 0.001,  pair_sig{k} = '***';
    elseif pair_pval(k) <= 0.01,   pair_sig{k} = '**';
    elseif pair_pval(k) <= 0.05,   pair_sig{k} = '*';
    else,                          pair_sig{k} = 'ns';
    end
end

% n/count labels above each bar
for i = 1:numBars
    rev_count = round(revs(i) * n_values(i));
    text(i, revs(i) + err_upper(i) + 0.02, sprintf('%d/%d', rev_count, n_values(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 8);
end

% Draw significance brackets, sorted by span (short brackets lower)
needed_ylim = 1.05;

spans = pair_list(:,2) - pair_list(:,1);
[~, sortIdx] = sort(spans);
pair_list_sorted = pair_list(sortIdx, :);
pair_sig_sorted  = pair_sig(sortIdx);
pair_pval_sorted = pair_pval(sortIdx);

bracketStep = 0.07;
baseHeight  = max(revs(:) + err_upper(:)) + 0.08;
tickLen     = 0.015;

% Greedy level assignment: track highest occupied level per bar position
levels = zeros(1, numBars);
bracketY = zeros(num_pairs, 1);
for k = 1:num_pairs
    x1 = pair_list_sorted(k, 1);
    x2 = pair_list_sorted(k, 2);
    lv = max(levels(x1:x2)) + 1;
    levels(x1:x2) = lv;
    bracketY(k) = baseHeight + (lv - 1) * bracketStep;
end

for k = 1:num_pairs
    x1 = pair_list_sorted(k, 1);
    x2 = pair_list_sorted(k, 2);
    y  = bracketY(k);
    plot([x1 x2], [y y],           'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
    plot([x1 x1], [y-tickLen y],   'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
    plot([x2 x2], [y-tickLen y],   'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
    text((x1+x2)/2, y+0.01, pair_sig_sorted{k}, ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
end
needed_ylim = max(bracketY) + 0.08;

% Print p-value table
fprintf('\n%s\n', repmat('=', 100, 1));
fprintf('Specified pairwise p-values (Fisher''s exact, two-sided, Holm-Bonferroni, k=%d):\n', num_pairs);
fprintf('%s\n', repmat('=', 100, 1));
fprintf('%-35s %-35s %10s %12s %10s\n', 'Group A', 'Group B', 'p-raw', 'p-adj (Holm)', 'Sig');
fprintf('%s\n', repmat('-', 100, 1));
for k = 1:num_pairs
    rp = raw_pval(k);
    ap = pair_pval(k);
    if rp < 0.001, rp_str = sprintf('%.2e', rp); else, rp_str = sprintf('%.4f', rp); end
    if ap < 0.001, ap_str = sprintf('%.2e', ap); else, ap_str = sprintf('%.4f', ap); end
    fprintf('%-35s %-35s %10s %12s %10s\n', conditions{pair_list(k,1)}, conditions{pair_list(k,2)}, rp_str, ap_str, pair_sig{k});
end
fprintf('%s\n', repmat('=', 100, 1));
fprintf('Significance levels: ns p>0.05, * p<=0.05, ** p<=0.01, *** p<=0.001, **** p<=0.0001\n');

% Formatting
xticks(1:numBars);
xticklabels(conditions);
xtickangle(45);
title('Longitudinal Muscle Ablation');
ylabel('Fraction inverted');
ylim([0, needed_ylim]);
set(gca, 'YTick', [0, 0.25, 0.5, 0.75, 1.0]);
set(gca, 'Box', 'off', 'TickDir', 'in', 'XColor', 'k');
set(gcf, 'color', 'w');
hold off;
legend('inverted', 'abnormal', 'location', 'northwest');
legend boxoff

saveas(gcf, './figure_longitudinalMuscleAblation.fig')
exportgraphics(gcf, './figure_longitudinalMuscleAblation.png', 'Resolution', 400)
exportgraphics(gcf, './figure_longitudinalMuscleAblation.pdf', 'ContentType', 'vector')
