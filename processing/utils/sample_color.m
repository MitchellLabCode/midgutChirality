function c = sample_color(idx, n_samples)
% sample_color  Color of sample idx in the per-embryo plots: WT in blues,
% byn>Myo1C in reds, from a coolwarm map with a gap between the two.
cpad = 3;
cmap = coolwarm(n_samples + cpad);
if mod(idx, 2)          % byn>Myo1C (odd)
    c = cmap(n_samples/2 + (idx-1)/2 + cpad + 1, :);
else                    % WT (even)
    c = cmap(idx/2 + 1, :);
end
end
