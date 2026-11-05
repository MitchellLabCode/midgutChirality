function [kappa, theta] = logrot(M)
    %%% logrot
    % Matrix logarithm of the 3D rotation matrix M
    % Inputs: M 3x3 array
    % Outputs: - kappa 3-vector (rotation axis, has norm theta)
    %          - theta float : rotation angle
    % NOTE that kappa is NOT a unit vector
    
    U = 0.5 * (M-M');
    kappa = [U(3,2), -U(3,1), U(2,1)];
    if norm(kappa) == 0
        kappa = [0., 0., 0.];
        theta = 0.;
    else
        nu = kappa/ norm(kappa);
        Kn = skewVec2Mat(nu);
        % for M = expm(theta*Kn): trace(Kn*M) = -2 sin(theta) and
        % trace(M) = 1 + 2 cos(theta). (Before 2026-10-01 the sine term was
        % missing its factor 1/2, giving atan(2 tan(theta)) ~ 2*theta.)
        theta = atan2(-0.5*trace(Kn*M), 0.5*(trace(M)-1));
        kappa = theta * nu;
    end
end

function U = skewVec2Mat(u)
    % Converts a 3x1 vector to a skew-symmetric matrix
    U = [  0     -u(3)   u(2);
          u(3)     0    -u(1);
         -u(2)    u(1)     0  ];
end