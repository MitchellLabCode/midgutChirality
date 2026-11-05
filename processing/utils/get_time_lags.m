function lags = get_time_lags(P, all_data)
% get_time_lags  Per-sample time lag [frames] used to align all embryos.
%
%   lags = get_time_lags(P, all_data)
%
% Reads the shifts written by timeAlign/timelag_fromRadii.m: per genotype,
% the shift of each embryo relative to its genotype's reference embryo, plus
% one cross-genotype shift for byn>Myo1C. The result is offset so that t = 0
% is the onset of the 2nd fold of the reference WT embryo (sample 4):
%
%   lags(i) = relative_shift(i) [+ cross_genotype_shift if byn>Myo1C] - fold_onset_2(sample 4)
%
% Aligned times are then built with aligned_time.m.

ts = load(P.timeshifts);
n = numel(P.samples);
lags = zeros(1, n);
for i = 1:n
    if P.samples(i).is_WT
        k = sum([P.samples(1:i).is_WT]);              % WT #k
        v = ts.relative_shifts.wt{k};
    else
        k = sum(~[P.samples(1:i).is_WT]);             % byn>Myo1C #k
        v = ts.relative_shifts.bynMyo1c{k};
    end
    if isempty(v), v = 0; end                         % the reference embryo
    if ~P.samples(i).is_WT, v = v + ts.cross_genotype_shift; end
    lags(i) = v;
end
fo = all_data{4}.fold_onset;
lags = lags - fo(2);
end
