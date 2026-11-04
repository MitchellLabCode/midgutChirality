% Fig 3a: Midgut Inversion Chirality Myo1C Overexpression Screen
% Plot results from Myo1COEScreen_Fig3a.xlsx
% Statistical test: two-sided Fisher's exact test with two-family Holm-Bonferroni
% Family 1 (midgut drivers, k=9) and Family 2 (hindgut drivers, k=7) corrected separately
% maternal nanos (condition 11) is excluded from both families
close all
clear all
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'))

% Read data from the Excel file
inputFile = fullfile(fileparts(mfilename('fullpath')), 'utils', 'Myo1COEScreen_Fig3a.xlsx');
data = readtable(inputFile, 'Sheet', 1, 'VariableNamingRule', 'preserve');

% Extract relevant data
conditions = data.Properties.VariableNames(2:end);
wt  = data{strcmp(data.("Midgut Coiling"), 'Percent WT'), 2:end};
rev = data{strcmp(data.("Midgut Coiling"), 'Percent inverted'), 2:end};
abnormal = 1 - wt - rev;
n_values = data{strcmp(data.("Midgut Coiling"), 'n'), 2:end};

% Reorder the conditions
order = [1, 3, 18, 5, 9, 17, 6, 7, 8, 11, 10, 12, 15, 2, 16, 4, 13, 14];
conditions = conditions(order);
wt       = wt(order);
rev      = rev(order);
abnormal = abnormal(order);
n_values = n_values(order);
numBars  = length(conditions);

% Create figure
figure; hold on;
b = bar([rev; abnormal]', 'stacked', 'EdgeColor', 'none');
b(1).FaceColor = [204/255, 51/255, 51/255];
b(2).FaceColor = [0.8, 0.8, 0.8];

% Wilson 95% CI error bars on inverted proportion
err_lower = zeros(numBars, 1);
err_upper = zeros(numBars, 1);
for i = 1:numBars
    rev_count = round(rev(i) * n_values(i));
    abn_count = round(abnormal(i) * n_values(i));
    wt_count  = round(wt(i)  * n_values(i));
    if (rev_count + abn_count + wt_count) > 0
        stats = threeway_uncertainty_wilson(rev_count, abn_count, wt_count);
        err_lower(i) = rev(i) - stats.ci95(1, 1);
        err_upper(i) = stats.ci95(1, 2) - rev(i);
    end
end
errorbar(1:numBars, rev, err_lower, err_upper, 'k', 'LineStyle', 'none', ...
    'CapSize', 5, 'LineWidth', 1.5, 'HandleVisibility', 'off');

%% Fisher's exact tests (two-sided) + two-family Holm-Bonferroni
% Family 1 — midgut-expressing drivers (conditions 2-10, k=9)
% maternal nanos (condition 11) — not in either family, uncorrected
% Family 2 — hindgut-expressing drivers (conditions 12-end, k=7)

n0  = n_values(1);
x0  = round(rev(1) * n0);

% Raw two-sided Fisher's exact p-values
raw_pvals = ones(numBars, 1);
for i = 2:numBars
    xi  = round(rev(i) * n_values(i));
    tbl = [xi, n_values(i)-xi; x0, n0-x0];
    [~, raw_pvals(i)] = fishertest(tbl, 'Tail', 'both');
end

% Two-family Holm-Bonferroni correction
fam1_idx = 2:10;
fam2_idx = 12:numBars;
k1 = length(fam1_idx);
k2 = length(fam2_idx);

pvals = ones(numBars, 1);

[sorted_p1, sort_idx1] = sort(raw_pvals(fam1_idx));
running_max1 = 0;
adj1 = zeros(k1, 1);
for r = 1:k1
    corrected1   = min(1, sorted_p1(r) * (k1 - r + 1));
    running_max1 = max(running_max1, corrected1);
    adj1(r)      = running_max1;
end
pvals(fam1_idx(sort_idx1)) = adj1;

[sorted_p2, sort_idx2] = sort(raw_pvals(fam2_idx));
running_max2 = 0;
adj2 = zeros(k2, 1);
for r = 1:k2
    corrected2   = min(1, sorted_p2(r) * (k2 - r + 1));
    running_max2 = max(running_max2, corrected2);
    adj2(r)      = running_max2;
end
pvals(fam2_idx(sort_idx2)) = adj2;

% maternal nanos (condition 11) — uncorrected. Neither a midgut nor hindgut driver, so not in either family. Report raw p-value only.
pvals(11) = raw_pvals(11);

sig_strs = repmat({''}, numBars, 1);
for i = 2:numBars
    if     pvals(i) <= 0.0001,  sig_strs{i} = '****';
    elseif pvals(i) <= 0.001,   sig_strs{i} = '***';
    elseif pvals(i) <= 0.01,    sig_strs{i} = '**';
    elseif pvals(i) <= 0.05,    sig_strs{i} = '*';
    end
end

% Labels above bars
for i = 1:numBars
    rev_count = round(rev(i) * n_values(i));
    if isempty(sig_strs{i})
        label_str = sprintf('\n%d/%d', rev_count, n_values(i));
    else
        label_str = sprintf('%s\n%d/%d', sig_strs{i}, rev_count, n_values(i));
    end
    text(i, rev(i) + err_upper(i) + 0.05, label_str, ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
end

% Report p-value table
fprintf('\n%s\n', repmat('=', 90, 1));
fprintf('Two-sided Fisher''s exact + two-family Holm-Bonferroni vs %s\n', conditions{1});
fprintf('Family 1 (midgut, k=%d): conditions 2-%d  |  Family 2 (hindgut, k=%d): conditions 12-%d  |  condition 11: uncorrected\n', k1, fam1_idx(end), k2, numBars);
fprintf('%s\n', repmat('=', 90, 1));
fprintf('%-40s %6s %12s %12s %10s\n', 'Condition', 'Family', 'p-raw', 'p-adj (Holm)', 'Sig');
fprintf('%s\n', repmat('-', 90, 1));
fprintf('%-40s %6s %12s %12s %10s\n', conditions{1}, '—', '—', '(reference)', '');
for i = 2:numBars
    if ismember(i, fam1_idx)
        fam_str = 'F1';
    elseif ismember(i, fam2_idx)
        fam_str = 'F2';
    else
        fam_str = '--';
    end
    if raw_pvals(i) < 0.001, rp = sprintf('%.2e', raw_pvals(i)); else, rp = sprintf('%.4f', raw_pvals(i)); end
    if pvals(i)     < 0.001, ap = sprintf('%.2e', pvals(i));     else, ap = sprintf('%.4f', pvals(i));     end
    fprintf('%-40s %6s %12s %12s %10s\n', conditions{i}, fam_str, rp, ap, sig_strs{i});
end
fprintf('%s\n', repmat('=', 90, 1));
fprintf('Significance: * p<=0.05, ** p<=0.01, *** p<=0.001, **** p<=0.0001\n');

% Formatting
xticks(1:numBars);
xticklabels(conditions);
xtickangle(45);
title('Midgut Inversion');
set(gca, 'YTick', [0, 0.25, 0.5, 0.75, 1.0]);
ylim([0, 1.05]);
ylabel('Fraction inverted');
set(gca, 'Box', 'off', 'TickDir', 'in', 'XColor', 'k');
set(gcf, 'color', 'w');
hold off;
legend('inverted', 'abnormal', 'location', 'northwest');
legend boxoff
saveas(gcf, './figure_Fig3a_OEScreen.fig')
exportgraphics(gcf, './figure_Fig3a_OEScreen.png', 'Resolution', 400)
exportgraphics(gcf, './figure_Fig3a_OEScreen.pdf', 'ContentType', 'vector')
