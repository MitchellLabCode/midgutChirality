function arrow3d(p, v, r_shaft, r_head, l_head, col)
% arrow3d  Solid 3D arrow from p along v: cylinder shaft + cone head, drawn
% as lit surfaces into the current axes (unlike quiver3's flat heads).
%
%   arrow3d(p, v, r_shaft, r_head, l_head, col)
L = norm(v); u = v / L;
l_head = min(l_head, 0.5*L);
n = 24; th = linspace(0, 2*pi, n+1);
tmp = [1 0 0]; if abs(u*tmp') > 0.9, tmp = [0 1 0]; end
e1 = cross(u, tmp); e1 = e1 / norm(e1);
e2 = cross(u, e1);
ring = @(r) r * (cos(th)' * e1 + sin(th)' * e2);
pts = @(R, z) R + z * u + p;
r0 = ring(r_shaft);
A = pts(r0, 0); B = pts(r0, L - l_head);
surf([A(:,1) B(:,1)], [A(:,2) B(:,2)], [A(:,3) B(:,3)], 'FaceColor', col, 'EdgeColor', 'none')
C = pts(ring(r_head), L - l_head);
tip = p + v;
surf([C(:,1) repmat(tip(1), n+1, 1)], [C(:,2) repmat(tip(2), n+1, 1)], ...
     [C(:,3) repmat(tip(3), n+1, 1)], 'FaceColor', col, 'EdgeColor', 'none')
c0 = p + (L - l_head) * u;
surf([C(:,1) repmat(c0(1), n+1, 1)], [C(:,2) repmat(c0(2), n+1, 1)], ...
     [C(:,3) repmat(c0(3), n+1, 1)], 'FaceColor', col, 'EdgeColor', 'none')
end
