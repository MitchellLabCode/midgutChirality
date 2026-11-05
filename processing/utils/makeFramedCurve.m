function [cent, Qs, kappa, omega] = makeFramedCurve(x, myocytes, meshes, n_cut, w, N_s)
%%% makeFramedCurve
% Builds the material frame of the midgut (Sec. IV.3 of the Supplement) and
% its Darboux vector kappa and spin omega.
%
% For every centerline point, the four myocyte stripes are intersected with
% the mesh cross-section at the same arclength (the two bracketing mesh hoops
% are blended), the intercept curves are smoothed in space and time, and d1 is
% the tangent-projected, inverse-square-weighted sum of the centerline-to-
% intercept directions over a window of +-n_ave points, with charges -1, -1,
% +1, +1 for the VL, DL, DR, VR stripes. d2 = d3 x d1, d3 = tangent.
% kappa (Eq. S14) and omega (Eq. S13) are matrix logarithms of the relative
% rotation between neighbouring frames (logrot.m); kappa is divided by the
% Voronoi length, omega is per frame.
%
% Makes the anatomical frame from the gut centerline and mycoyte curves.
% Inputs:
%  - x : n_times x (n_points+n_cut) x 3 array containing the centerline
%  - myocytes: struct with myocyte curves (n_times x n_points x 3 arrays).
%   - mesh : n_timesx1  struct containing vertices and faces information
%   with nU = 100, nV = 100 vertices
%  - n_cut: int specific how many points to ignore at the beginning of
%  centerline
% Outputs:
% cent : (n_times)x(n_points)x3 centerline cropped by n_cut, and interpolated to match
%  Q : (n_times)x(n_points-1)x3x3 array representing the material frame
% %note - Q(t, i, 1,:) is m1, Q(t, i,2,:) is m2, Q(t, i,3,:) is the tangent
% % at time t and position i.
% kappa : n_times x(n_points-2)x3  Darboux vector in the lab frame
% omega : (n_times-1)x(n_points-1)x3 Rotation vector in the lab frame
% the myocyte curves

%%% definition of material frames
% tangent : (n_points-1)x3 material tangent - points A to P
% m1_bar : (n_points-1)x3 material frame m1 - points L to R 
% m2_bar : (n_points-1)x3 material frame m2 - points V to D ('up')

intercept_fail = false;
[n_times, n_points, ~] = size(myocytes.red);

% prepare the arrays
cent = zeros(n_times, n_points,3);
Qs = zeros(n_times, n_points-1,3,3);
kappa = zeros(n_times, n_points-2,3);
omega = zeros(n_times-1, n_points-1,3);

curves_intercepts= zeros(n_times, n_points,3, 4);
ell = 0.3;
tau = 0.3;

coeff_j = [-1,-1,1,1]; % red, green, yellow, blue

for t = 1:n_times

    curves = {myocytes.red(t,:,:), myocytes.green(t,:,:), myocytes.yellow(t,:,:), myocytes.blue(t,:,:)};

    msh  = meshes{t};

    vert_UV = reshape(msh.v, [100,100,3]);
    % centerline interpolation - cut the beginning
    %n_cent = size(x,2);

    cent_c = squeeze(x(t,n_cut+1:end, :));

    tangent = cent_c(2:end,:) - cent_c(1:end-1,:);
    tangent = tangent ./ vecnorm(tangent,2,2);
    
    %curves_intercepts = zeros(size(cent_c,1), 3, 4);
    % for each myocyte curve, compute the associated ribbon frame

    %m1_bar = zeros(size(tangent,1), 3);

    % --- V6 FIX: arclength correspondence centerline point -> mesh hoop ----
    % cent_c is x(t, n_cut+1:end, :), and x is resampled uniformly in
    % arclength over the full curve, while the mesh hoops are spaced
    % non-uniformly. So centerline point k does not live on hoop k. Build
    % the mapping by matching arclength along the hoop-centroid curve
    % (which is what the raw TubULAR centerline avgpts actually are).
    hoop_c  = squeeze(mean(vert_UV, 2));                        % (n_hoop,3) centroids
    n_hoop  = size(hoop_c, 1);
    hoop_ss = [0; cumsum(vecnorm(diff(hoop_c, 1, 1), 2, 2))];   % (n_hoop,1)

    x_full = squeeze(x(t,:,:));                                 % uncut centerline
    x_ss   = [0; cumsum(vecnorm(diff(x_full, 1, 1), 2, 2))];
    s_cent = x_ss(n_cut+1:end);                                 % arclength of cent_c pts
    if x_ss(end) > 0
        s_cent = s_cent * (hoop_ss(end) / x_ss(end));           % smoothing shortens x
    end

    % fractional hoop index for each centerline point
    hoop_of_k = interp1(hoop_ss, (1:n_hoop)', s_cent, 'linear', 'extrap');
    hoop_of_k = min(max(hoop_of_k, 1), n_hoop);

    for j = 1:4
        c_curv = curves{j};
        for k = 1:(size(cent_c,1)-1)
            % blend the two hoops bracketing this point's arclength, so the
            % intercept is taken at the same location as cent_c(k,:)
            mf  = hoop_of_k(k);
            mi0 = floor(mf);
            mi1 = min(mi0 + 1, n_hoop);
            aa  = mf - mi0;
            ring = (1-aa) * squeeze(vert_UV(mi0,:,:)) ...
                 +    aa  * squeeze(vert_UV(mi1,:,:));
            intercept_point = UVplaneintercept(ring, squeeze(c_curv));
            if ~isempty(intercept_point)
                curves_intercepts(t,k, :, j) = intercept_point;
                %m1 =  intercept_point - cent_c(k,:);
                %m1 = m1 ./ norm(m1,2);
                %m1 = m1 - dot(m1, squeeze(tangent(k,:))) * squeeze(tangent(k,:));
                %m1_bar(k,:) = m1_bar(k,:) + 0.5*coeff_j(j) * m1;

            else
                intercept_fail = true;
                curves_intercepts(t,k, :, j) = cent_c(k,:);
            end
        end
        end


    % s_interpolation? 
    %cent_c = [c_x; c_y; c_z]';
    cent(t,:,:) = cent_c;
    %m1_bar = m1_bar ./ vecnorm(m1_bar, 2, 2);
    %m2_bar = cross(tangent, m1_bar, 2);

   

   % make the adapted anatomical frame 'LR-DV-AP' - y, z, x is AP, 
   % Q(3,:) is along the tangent, Q(1,:) is along LR (so encodes dorso-vental
   % rotation)
   % and Q(2,:) is along DV (so encodes LR rotations)

   if isnan(tangent(1,:))
       disp('Nan!')
   end

   %Qs(t,:,1,:) = m1_bar;
   %Qs(t,:,2,:) = m2_bar;
   %Qs(t,:,3,:) = tangent;

   %lengths_i = vecnorm( cent_c(2:end,:)-cent_c(1:end-1,:),2,2);


   %for i = 1:size(Qs,1)-1
   %     Q1 = squeeze(Qs(t,i,:,:));
   %     Q2 = squeeze(Qs(t,i+1,:,:));
   %     [kappa_i, ~] = logrot(Q1' * Q2);
   %     voro = 0.5*(lengths_i(i)+lengths_i(i+1));
   %     kappa(t,i,:) = kappa_i / voro;
        %dar_vec(i,1) = dot(kappa_i, m1_bar(i,:), 2);
        %dar_vec(i,2) = dot(kappa_i, m2_bar(i,:), 2);
        %dar_vec(i,3) = dot(kappa_i, tangent(i,:), 2);
   %end
end

% now smooth curves
equi = false;
for j=1:4
curves_intercepts(:,:,:,j) = smooth_stdata(squeeze(curves_intercepts(:,:,:,j)), ell, tau, N_s, equi);
end

n_ave = 5;
% now compute m1, m2, and 
for t= 1:n_times
    tangent = squeeze(cent(t,2:end,:) - cent(t,1:end-1,:));
    tangent = tangent ./ vecnorm(tangent,2,2);
    
    %curves_intercepts = zeros(size(cent_c,1), 3, 4);
    % for each myocyte curve, compute the associated ribbon frame
    
    m1_bar = zeros(size(tangent,1), 3);
    num_q = size(tangent,1);
    for j =1:4
        
        for k = 1:size(tangent,1)
            m1 = zeros(3,1);
            % average 'electric field'

            left_idx = max(1,k-n_ave);
            right_idx =  min(num_q,k+n_ave);
            ci_c = squeeze(curves_intercepts(t,left_idx:right_idx,:,j));
            m0 = ci_c - squeeze(cent(t,left_idx:right_idx,:));
            is_valid = vecnorm(m0,2,2) > 1e-2;
            m0 = m0(is_valid, :);
            m1 = sum(m0 ./ (vecnorm(m0,2,2) .^3),1);

            
            m1 = m1 - dot(m1, squeeze(tangent(k,:)),2) * squeeze(tangent(k,:));
            %m1 = m1 ./ vecnorm(m1,2,2);
            m1_bar(k,:) = m1_bar(k,:) + 0.5* coeff_j(j) * m1;
        end
        
    end
    m1_bar = m1_bar ./ vecnorm(m1_bar, 2, 2);
    m2_bar = cross(tangent, m1_bar, 2);

   Qs(t,:,1,:) = m1_bar;
   Qs(t,:,2,:) = m2_bar;
   Qs(t,:,3,:) = tangent;

   % segment lengths of THIS timepoint's centerline (cent_c is left over
   % from the loop above and holds the last timepoint only)
   lengths_i = vecnorm( squeeze(cent(t,2:end,:)-cent(t,1:end-1,:)),2,2);


   % one kappa per interior vertex: size(Qs,2) is the number of edges
   for i = 1:size(Qs,2)-1
        Q1 = squeeze(Qs(t,i,:,:));
        Q2 = squeeze(Qs(t,i+1,:,:));
        [kappa_i, ~] = logrot(Q1' * Q2);
        voro = 0.5*(lengths_i(i)+lengths_i(i+1));
        kappa(t,i,:) = kappa_i / voro;
        %dar_vec(i,1) = dot(kappa_i, m1_bar(i,:), 2);
        %dar_vec(i,2) = dot(kappa_i, m2_bar(i,:), 2);
        %dar_vec(i,3) = dot(kappa_i, tangent(i,:), 2);
   end
end
    

% compute omega
for t = 1:(n_times-1)
    for i = 1:size(Qs,2)
        Q1 = squeeze(Qs(t,i,:,:));
        Q2 = squeeze(Qs(t+1,i,:,:));
        [om, ~] = logrot(Q1' * Q2);
        omega(t,i,:) = om;
    end
end
if intercept_fail
    warning('Interception failure. Increase n_cut')
end

end