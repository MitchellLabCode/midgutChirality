function xs = smooth_stdata(x, ell, tau, N, equi)
% Smoothes a 3D space curve in space and time. Uses Taubin smoothing in
% space, gaussian smoothing in time.
% returns an arclength parameterized curve
% Inputs
% x : num_t x num_s x 3 array
% ell: Taubin parameter
% tau: low pass time scale
% N: number of Taubin iterations
% Outputs
% xs : num_t x num_s x 3 array - smoothed curve.
%
%  Taubin, Curve and surface smoothing without shrinkage (1995)

if nargin == 4
    equi = true;
end

[num_t, num_s, ~] = size(x);

% Reparameterize to arclength parameter and
% compute lambda_x(t) mu_x(t), lambda_t, mu_t.

if equi
lambda_x = zeros(num_t, 1);
mu_x = zeros(num_t, 1);
xs = zeros(size(x));
for t = 1:num_t
    x_t  = squeeze(x(t,:,:));
    ds = vecnorm(x_t(2:end,:)-x_t(1:end-1,:), 2,2);
    s_arr = [0.; cumsum(ds)];
    L = s_arr(end);
    s_interp = linspace(0, L, num_s);
    if t == 1
        L0 = L;
    end

    lambda_x(t) = ell;%*(L0/L)^2;
    mu_x(t) = -(ell+2e-2);%*(L0/L)^2;

    %[S_u,idx_u,~] = unique(s_arr);
    xs(t, :, 1) = interp1(s_arr, x_t(:,1), s_interp);
    xs(t, :, 2) = interp1(s_arr, x_t(:,2), s_interp);
    xs(t, :, 3) = interp1(s_arr, x_t(:,3), s_interp);

    %xs(t, :, 1) = interp1(S_u, x_t(idx_u,1), s_interp);
    %xs(t, :, 2) = interp1(S_u, x_t(idx_u,2), s_interp);
    %xs(t, :, 3) = interp1(S_u, x_t(idx_u,3), s_interp);

end
else
    xs = x; % no interpolation
    lambda_x = ell*ones(num_t, 1);
    mu_x = -(ell+2e-2)*ones(num_t, 1);
end


% apply Taubin smoothing in space but Gaussian smoothing in time.



for t = 1:num_t

    ell_t = lambda_x(t);
    mu_t = mu_x(t);
    for n = 1:N
        xs(t, 2:end-1, :) = xs(t, 2:end-1, :) + ell_t *(xs(t, 1:end-2,:) + xs(t, 3:end,:) - 2*xs(t, 2:end-1, :));
        xs(t, 2:end-1, :) = xs(t, 2:end-1, :) + mu_t *(xs(t, 1:end-2,:) + xs(t, 3:end,:) - 2*xs(t, 2:end-1, :));
    end

end


% time smoothing
ell_t = tau;
mu_t = 0*-(tau-1e-2);
for i = 1:num_s
    %xs(:,i,:) = smoothdata(squeeze(x(:,i,:)), 2, 'sgolay', tau);

    for n = 1:N
        xs(2:end-1, i, :) = xs(2:end-1,i, :) + ell_t *(xs(1:end-2,i, :) + xs(3:end,i,:) - 2*xs(2:end-1,i, :));
        %xs(2:end-1, i, :) = xs(2:end-1,i, :) + mu_t *(xs(1:end-2,i,:) + xs(3:end,i, :) - 2*xs(2:end-1,i, :));
    end
end

end