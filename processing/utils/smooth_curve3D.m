function smoothed_curve = smooth_curve3D(init_curve, w, n_c, tol, max_iters)
%smoothed_curve = smooth_curve3D(init_curve, w, n_c, tol, max_iters)
%   Diffusion-like flow of a curve towards a target geometry.
%
% Parameters 
% ----------
% init_curve : Nx1 3D float
%   initial curve
% w : float 
%   weight for curvature-based flow - bigger means more smoothing relative
%   to staying 
% n_c : int
%    #points along the curve for resampling
% tol : float
%   tolerance for convergence of the flow: fraction of E=(summed distance
%   to the closest matching points + w*overall length of the curve)
% max_iters : int
%   maximum # iterations to permit for convergence
%
% Requirements 
% ------------
% distance2curve : function in matlab/curve_handling
%
% Example Usage 
% -------------
%
%
%
%
% Nico Romeo 2025, inspired by Santana-Cedrés et al, IPOL 2021.,
% https://www.ipol.im/pub/art/2021/292/?utm_source=doi


% get arc-length of the input curve - rescale all coordinates by initial
% arc-length.

L_init = sum(vecnorm(init_curve(2:end,:)-init_curve(1:end-1,:),2,2));
ic_scaled = init_curve / L_init;

% equi-distant sampling
n_c = n_c - 1 ;
uni_curve = resampleCurve(ic_scaled, 1/n_c);

% 
delta = min([1/(2*n_c), (1/n_c^2)/(pi*w)]);
cur_curve = [uni_curve; ic_scaled(end,:)];


% Check that the curve is the right size
if size(cur_curve, 1) ~= n_c + 1
    % Repeat the first row of cur_curve to match the required size
    cur_curve = [repmat(cur_curve(1, :), n_c + 1 - size(cur_curve, 1), 1); cur_curve];
end
   
assert(size(cur_curve, 1) == n_c + 1) 

% Check that no imaginary component has appeared
%assert(all(imag(cur_curve)<1e-16))

% until convergence or maxiters:
Em1 =0.0;
for i=1:max_iters

    try
        % Check that no imaginary component has appeared
        assert(all(sum(imag(cur_curve))<1e-10))
    catch
        error(['imaginary component appeared on iteration=' num2str(i)])
    end


    %En = 0.0; % energy measure
    % uni_curve=real(uni_curve);
    % cur_curve=real(cur_curve);

    [closest_points, dist_closest,~] = distance2curve(uni_curve, cur_curve);

    normals = 0.5*(cur_curve(1:end-2,:) + cur_curve(3:end,:)) - cur_curve(2:end-1,:);
    norms_normals = vecnorm(normals, 2, 2);
    normals = normals ./ norms_normals;
    normals(norms_normals < 1e-5, :) = 0; 
    new_curve = zeros(size(cur_curve));
    new_curve(1,:) = cur_curve(1,:);
    new_curve(end,:) = cur_curve(end,:);

    for n = 2:size(cur_curve,1)-1
        % for each interior point : compute distance to curve and point of intersection,
        % curve shortening flow
        % update positions
        % update energy measure
        tan1 = (cur_curve(n+1,:)- cur_curve(n,:)) / norm(cur_curve(n+1,:)- cur_curve(n,:));
        tan2 = (cur_curve(n,:)- cur_curve(n-1,:)) / norm(cur_curve(n,:)- cur_curve(n-1,:));
    
        % curvature of the curve at this location is arccos(tan1.tan2) * ds
        kappa = acos(dot(tan1, tan2))*n_c ;  
        if isnan(kappa)
            kappa = 0.0;
        end
        if dist_closest(n) > 1e-9
            new_curve(n,:) = cur_curve(n,:) - delta * (cur_curve(n,:) - closest_points(n,:)) /dist_closest(n) +  delta* w* kappa * normals(n-1,:);
        else
            new_curve(n,:) = cur_curve(n,:) +  delta*w* kappa * normals(n-1,:);
        end
    end
    % resample
    %cur_curve = resampleCurve(new_curve, 1/n_c);
    %cur_curve = [cur_curve; new_curve(end,:)];
    cur_curve = new_curve;
    En = sum(dist_closest) + w*size(cur_curve,1);
    if abs(Em1-En)/(Em1+1e-16) < tol
        break
    end

    Em1 = En;

    end
    smoothed_curve = L_init * cur_curve;
end



