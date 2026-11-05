function D = compute_deformation_measures(d, lag, P)
% compute_deformation_measures  Fig. 2 deformation measures of one embryo
% (Sec. IV.4 of the Supplement), per region and time.
%
%   D = compute_deformation_measures(d, lag, P)
%
% d   : one element of all_data (all_framesInterpV6_10samples.mat)
% lag : time lag of this sample [frames] (get_time_lags.m)
%
% Every measure is an (n x 3) array, columns = [anterior, posterior, all_midgut]:
%   D.E.t                          aligned time [min] (frames)
%   D.E.E_b1, E_b2, E_b, E_tw      mean of (k1^2)/2, (k2^2)/2, (k1^2+k2^2)/2, (k3^2)/2   [um^-2]
%   D.E.E_b1_total ... E_tw_total  the same integrated over the region                  [um^-1]
%   D.E.Theta_b1, _b2, _b, _tw     total angle: integral of |k1|, |k2|, |(k1,k2)|, |k3|  [rad]
%   D.omega.t                      aligned time [min] (intervals)
%   D.omega.omega_b, _tw, _1, _2   mean of |(w1,w2)|, |w3|, |w1|, |w2|                    [rad/min]
%   D.mu.t                         aligned time [min] (intervals)
%   D.mu.mu_b, _tw, _1, _2         mean of |(dk1,dk2)|, |dk3|, |dk1|, |dk2| (time deriv.) [rad/(um min)]
% with ki = kappa . d_i and wi = omega . d_i the material components.
%
% Numerics (identical for all measures): 3x3 median filter in the (time,
% arclength) plane; the first/last time row and arclength column, where the
% filter is undefined, are dropped; Voronoi-length weights for kappa and
% d_t kappa (vertices), edge-length weights for omega (edges); the
% anterior/posterior split is the time-dependent max-x point of the
% centerline. omega and d_t kappa are per frame divided by the frame interval.

dt = P.dt;
Q = d.Q; x = d.cent; kap = d.kappa; om = d.omega;
num_t = size(Q, 1);

ell  = vecnorm(x(:,2:end,:) - x(:,1:end-1,:), 2, 3);      % (num_t, 99) edge lengths
voro = 0.5 * (ell(:,1:end-1) + ell(:,2:end));              % (num_t, 98) Voronoi lengths
idx_split = zeros(num_t, 1);
for t = 1:num_t, [~, idx_split(t)] = max(x(t,:,1)); end

%%% curvature measures E and total angles Theta (on frames)
k1 = medfilt2(squeeze(dot(kap, squeeze(Q(:,1:end-1,1,:)), 3)), [3 3]);
k2 = medfilt2(squeeze(dot(kap, squeeze(Q(:,1:end-1,2,:)), 3)), [3 3]);
k3 = medfilt2(squeeze(dot(kap, squeeze(Q(:,1:end-1,3,:)), 3)), [3 3]);
tv = 2:(num_t-1); sv = 2:97;                               % kappa column s = point s+1
k1 = k1(tv,sv); k2 = k2(tv,sv); k3 = k3(tv,sv); w = voro(tv,sv);
mA = ap_mask(sv, idx_split(tv) - 1);
D.E.t = aligned_time(tv(:), lag, 'frame0', P);
[D.E.E_b1, D.E.E_b1_total] = region_stats(0.5*k1.^2, w, mA);
[D.E.E_b2, D.E.E_b2_total] = region_stats(0.5*k2.^2, w, mA);
[D.E.E_b,  D.E.E_b_total]  = region_stats(0.5*(k1.^2 + k2.^2), w, mA);
[D.E.E_tw, D.E.E_tw_total] = region_stats(0.5*k3.^2, w, mA);
[~, D.E.Theta_b1] = region_stats(abs(k1), w, mA);
[~, D.E.Theta_b2] = region_stats(abs(k2), w, mA);
[~, D.E.Theta_b]  = region_stats(sqrt(k1.^2 + k2.^2), w, mA);
[~, D.E.Theta_tw] = region_stats(abs(k3), w, mA);

%%% rotation rates omega (on intervals f -> f+1; frame taken at f+1)
w1 = medfilt2(squeeze(dot(om, squeeze(Q(2:end,:,1,:)), 3)), [3 3]);
w2 = medfilt2(squeeze(dot(om, squeeze(Q(2:end,:,2,:)), 3)), [3 3]);
w3 = medfilt2(squeeze(dot(om, squeeze(Q(2:end,:,3,:)), 3)), [3 3]);
wv = 2:(size(om,1) - 1); sv = 2:98;                        % edge e joins points e, e+1
w1 = w1(wv,sv); w2 = w2(wv,sv); w3 = w3(wv,sv); w = ell(wv,sv);
mA = ap_mask(sv, idx_split(wv));
D.omega.t = aligned_time(wv(:), lag, 'interval0', P);
D.omega.omega_b  = region_stats(sqrt(w1.^2 + w2.^2)/dt, w, mA);
D.omega.omega_tw = region_stats(abs(w3)/dt, w, mA);
D.omega.omega_1  = region_stats(abs(w1)/dt, w, mA);
D.omega.omega_2  = region_stats(abs(w2)/dt, w, mA);

%%% curvature change mu (on intervals f -> f+1)
c1 = squeeze(dot(kap, squeeze(Q(:,1:end-1,1,:)), 3));
c2 = squeeze(dot(kap, squeeze(Q(:,1:end-1,2,:)), 3));
c3 = squeeze(dot(kap, squeeze(Q(:,1:end-1,3,:)), 3));
g1 = medfilt2(c1(2:end,:) - c1(1:end-1,:), [3 3]);
g2 = medfilt2(c2(2:end,:) - c2(1:end-1,:), [3 3]);
g3 = medfilt2(c3(2:end,:) - c3(1:end-1,:), [3 3]);
dv = 2:(size(g1,1) - 1); sv = 2:97;
g1 = g1(dv,sv); g2 = g2(dv,sv); g3 = g3(dv,sv);
vm = 0.5 * (voro(1:end-1,:) + voro(2:end,:));              % Voronoi length averaged over f, f+1
w = vm(dv,sv);
mA = ap_mask(sv, idx_split(dv) - 1);
D.mu.t = aligned_time(dv(:), lag, 'interval0mid', P);
D.mu.mu_b  = region_stats(sqrt(g1.^2 + g2.^2)/dt, w, mA);
D.mu.mu_tw = region_stats(abs(g3)/dt, w, mA);
D.mu.mu_1  = region_stats(abs(g1)/dt, w, mA);
D.mu.mu_2  = region_stats(abs(g2)/dt, w, mA);
end


function mA = ap_mask(sv, split)
% anterior: columns before the max-x point (one row per kept time)
mA = false(numel(split), numel(sv));
for k = 1:numel(split), mA(k,:) = sv < split(k); end
end

function [mn, tot] = region_stats(F, w, mA)
% per row: weighted mean and integral over [anterior, posterior, all midgut]
n = size(F, 1); mn = zeros(n, 3); tot = zeros(n, 3);
for t = 1:n
    a = mA(t,:); p = ~a; ww = w(t,:); ff = F(t,:);
    tot(t,:) = [sum(ww(a).*ff(a)), sum(ww(p).*ff(p)), sum(ww.*ff)];
    L = [sum(ww(a)), sum(ww(p)), sum(ww)];
    mn(t,:) = tot(t,:) ./ L; mn(t, L == 0) = 0;
end
end
