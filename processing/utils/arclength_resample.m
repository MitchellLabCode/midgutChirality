function C = arclength_resample(C)
% arclength_resample  Resample every frame of an (n_t x N x 3) curve array to
% N points equally spaced in arclength.
for t = 1:size(C, 1)
    ce = squeeze(C(t,:,:));
    s = [0; cumsum(vecnorm(diff(ce), 2, 2))];
    C(t,:,:) = interp1(s / s(end), ce, linspace(0, 1, size(C, 2)), 'linear');
end
end
