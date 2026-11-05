function lbl = sample_label(P, i)
% sample_label  Display name of sample i: 'WT #k' or 'byn>Myo1C #k', with k
% the embryo's rank within its genotype (sample 4 -> 'WT #2').
if P.samples(i).is_WT
    lbl = sprintf('WT #%d', sum([P.samples(1:i).is_WT]));
else
    lbl = sprintf('byn>Myo1C #%d', sum(~[P.samples(1:i).is_WT]));
end
end
