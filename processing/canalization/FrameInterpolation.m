function Q = FrameInterpolation(xi, Q1, Q2, tau)
%FRAMEINTERPOLATION Interpolate adapted frames between two centerlines.
%   Q = FrameInterpolation(xi, Q1, Q2, tau)
%
% Frames are stored as rows: Q(s,i,:) is the i-th frame vector of segment s,
% and Q(s,3,:) is the unit tangent. A rotation R acting on the frame vectors
% is therefore applied as Qs * R.'.
%
% Inputs:
%   xi  : [N x 3] interpolated centerline
%   Q1  : [N-1 x 3 x 3] frames of curve 1 (tau = 0, WT)
%   Q2  : [N-1 x 3 x 3] frames of curve 2 (tau = 1, Myo)
%   tau : interpolation parameter
%
% Steps (Fig. S15b):
%   i.   minimally rotate Q2 so its tangent matches the tangent of Q1
%   ii.  interpolate about the shared tangent: rotate Q1 by tau*phi, where
%        phi is the angle about the tangent from Q1 to the rotated Q2
%   iii. minimally rotate the result onto the tangent of xi
% Every step is a rotation, so Q stays orthonormal and right-handed.

num_s = size(Q1,1);
Q = zeros(size(Q1));
for s = 1:num_s
    Q1a = squeeze(Q1(s,:,:));
    Q2a = squeeze(Q2(s,:,:));
    t1 = Q1a(3,:).';

    % i. bring Q2 into the adapted plane of Q1
    R1 = minimal_rotation(Q2a(3,:), t1);
    Q2a = Q2a * R1.';

    % ii. interpolate the angle about the shared tangent
    phi = atan2(dot(cross(Q1a(1,:), Q2a(1,:)), t1), dot(Q1a(1,:), Q2a(1,:)));
    Q3a = Q1a * axis_rotation(t1, tau * phi).';

    % iii. bring to the adapted plane of xi
    t3 = xi(s+1,:) - xi(s,:);
    t3 = t3 / norm(t3, 2);
    R2 = minimal_rotation(t1, t3);

    Q(s,:,:) = Q3a * R2.';
end

end

function R = axis_rotation(u, theta)
%AXIS_ROTATION Rotation by angle theta about unit axis u (Rodrigues).
u = u(:) / norm(u);
K = [   0   -u(3)  u(2);
      u(3)   0    -u(1);
     -u(2)  u(1)   0   ];
R = eye(3) + sin(theta) * K + (1 - cos(theta)) * (K * K);
end

function R = minimal_rotation(a, b)
%MINIMAL_ROTATION Minimal 3D rotation mapping unit vector a to unit vector b.
% R = minimal_rotation(a, b)
%
% Inputs:
% a, b : 3x1 or 1x3 vectors (expected unit; normalized internally)
%
% Output:
% R : 3x3 rotation matrix with R*a = b, using the shortest rotation.
%
% Handles edge cases:
% - a ~ b -> R = I
% - a ~ -b -> 180-degree rotation about any axis orthogonal to a

a = a(:); b = b(:);
if numel(a) ~= 3 || numel(b) ~= 3
    error('a and b must be 3D vectors.');
end

na = norm(a); nb = norm(b);
if na < eps || nb < eps
    error('a and b must be non-zero.');
end
a = a / na;
b = b / nb;

c = dot(a, b);                % cos(theta)
c = max(min(c, 1.0), -1.0);   % numerical clamp

% Parallel case
if c > 1 - 1e-12
    R = eye(3);
    return;
end

% Antiparallel case (180 deg): choose stable orthogonal axis
if c < -1 + 1e-12
    % pick basis vector least aligned with a
    [~, idx] = min(abs(a));
    e = zeros(3,1); e(idx) = 1;
    u = cross(a, e);
    u = u / norm(u);  % unit axis orthogonal to a

    % Rotation by pi: R = -I + 2*u*u^T
    R = -eye(3) + 2*(u*u.');
    return;
end

% General case: Rodrigues with v = a x b
v = cross(a, b);
s = norm(v);  % sin(theta)

K = [   0   -v(3)  v(2);
      v(3)   0    -v(1);
     -v(2)  v(1)   0   ];

% Equivalent to Rodrigues: R = I + K + K^2 * ((1-c)/s^2)
R = eye(3) + K + K*K * ((1 - c) / (s^2));

end
