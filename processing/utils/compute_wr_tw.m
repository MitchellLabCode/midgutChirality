function [Wr, Tw, L] = compute_wr_tw(x, Q)
% compute_wr_tw  Writhe, twist and length of a framed curve over time.
%
%   [Wr, Tw, L] = compute_wr_tw(x, Q)
%
% x : (n_t x N x 3) centerline;  Q : (n_t x N-1 x 3 x 3) frame (d1 = Q(:,:,1,:))
% Wr: writhe of the centerline (compute_writhe.m); Tw: twist of d1 about it
% (compute_twist.m); L: centerline length.

n_t = size(x, 1);
Wr = zeros(n_t, 1); Tw = zeros(n_t, 1); L = zeros(n_t, 1);
for t = 1:n_t
    xt = squeeze(x(t,:,:));
    L(t)  = sum(vecnorm(diff(xt, [], 1), 2, 2));
    Tw(t) = compute_twist(xt, squeeze(Q(t,:,1,:)));
    Wr(t) = compute_writhe(xt);
end
end
