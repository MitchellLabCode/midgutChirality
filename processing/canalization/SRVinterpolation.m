function z = SRVinterpolation(x, y, tau)
    %SRVinterpolation SRV geodesic interpolation between two discrete curves.
    %   z = SRVinterpolation(x, y, tau)
    %
    % Inputs:
    %   x, y : [N x D] sampled curves (same size)
    %   tau  : scalar or vector in [0,1]
    %
    % Output:
    %   z    : if tau scalar -> [N x D]
    %          if tau vector -> [N x D x numel(tau)]

    if ~isequal(size(x), size(y))
        error('x and y must have the same size.');
    end
    if size(x,1) < 2
        error('Curves must have at least 2 points.');
    end

    qx = xtoq_local(x);
    qy = xtoq_local(y);

    d = dot(qx(:), qy(:));
    d = min(max(d, -1.0), 1.0);
    theta = acos(d);

    if isscalar(tau)
        z = one_interp(qx, qy, theta, tau);
    else
        tau = tau(:).';
        nT = numel(tau);
        z = zeros(size(x,1), size(x,2), nT);
        for k = 1:nT
            z(:,:,k) = one_interp(qx, qy, theta, tau(k));
        end
    end
end

function z = one_interp(qx, qy, theta, tau)
    if theta < 1e-12
        alpha = qx;  % avoid 0/0 when curves are nearly identical in SRV space
    else
        alpha = (sin(theta*(1-tau)) * qx + sin(theta*tau) * qy) / sin(theta);
    end
    z = qtox_local(alpha, zeros(1, size(qx,2))); % matches your Python default
end

function q = xtoq_local(x)
% Python equivalent:
% tangent = x[1:,:] - x[:-1,:]
% norm_tangent = np.linalg.norm(tangent,axis=1)
% L = np.sum(norm_tangent)
% return tangent / np.sqrt(norm_tangent[:,None]*L)

    tangent = diff(x, 1, 1);                       % [N-1 x D]
    norm_tangent = sqrt(sum(tangent.^2, 2));       % [N-1 x 1]
    L = sum(norm_tangent);

    denom = sqrt(norm_tangent * L);                % [N-1 x 1]
    denom(denom == 0) = eps;                       % avoid divide-by-zero
    q = tangent ./ denom;                          % implicit expansion
end

function x = qtox_local(q, x0)
% Python equivalent:
% normq = np.linalg.norm(q, axis=-1)
% tangent = q * normq[:,None]
% x = np.vstack([zeros(3), cumsum(tangent,axis=0)])
% return x0 + x

    if nargin < 2
        x0 = zeros(1, size(q,2));
    end

    normq = sqrt(sum(q.^2, 2));                    % [N-1 x 1]
    tangent = q .* normq;                          % [N-1 x D]
    x = [zeros(1, size(q,2)); cumsum(tangent, 1)];
    x = x + x0;
end

function d = elasticmetric_srv(x, y)
% Optional helper matching your Python elasticmetric().
    qx = xtoq_local(x);
    qy = xtoq_local(y);
    c = dot(qx(:), qy(:));
    if c > 1.0
        warning('warning above 1: %.16g', c);
        d = 0.0;
    elseif c < -1.0
        warning('warning below -1: %.16g', c);
        d = pi;
    else
        d = acos(c);
    end
end