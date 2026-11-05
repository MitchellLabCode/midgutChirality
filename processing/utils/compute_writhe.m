function Wr = compute_writhe(centerline)
  % 
  %   This function computes the total writhe history of a rod.
  %   Equations used are from method 1a from Klenin & Langowski 2000 paper.
  % 
  %   Parameters
  %   ----------
  %   center_line : 
  %       3D (n_nodes, 3) array containing data with 'float' type.
  %       rod node positions.
  % 
  %   Returns
  %   -------
  %     Wr, float: writhe of curve
  %     
  [N, ~] = size(centerline);


  omega_star = zeros(N- 2, N - 1);
  segment_wr = zeros(N-2, N-1);

  for i = 1:(N-2)
      for j = (i+1):(N-1)
            point_one = centerline(i,:);
            point_two = centerline(i+1,:);
            point_three = centerline(j,:);
            point_four = centerline(j+1,:);

            r12 = point_two -point_one;
            r34 = point_three -point_four;
            r14 = point_four -point_one;
            r13 = point_three -point_one;
            r23 = point_three -point_two;
            r24 = point_four -point_two;

            n1 = cross(r13, r14);
            n1 = n1/ norm(n1);

            n2 = cross(r14, r24);
            n2 = n2/ norm(n2);

            n3 = cross(r24, r23);
            n3 = n3/ norm(n3);
            
            n4 = cross(r23, r12);
            n4 = n4/ norm(n4);

            omega_star(i,j) = asin(dot(n1,n2)) ...
                            + asin(dot(n2,n3)) ...
                            + asin(dot(n3,n4)) ...
                            + asin(dot(n4,n1)) ;
            if isnan(omega_star(i,j))
                omega_star(i,j) = 0.;
            end
            segment_wr(i,j) = -omega_star(i,j) * sign(dot(cross(r34,r12),r13))/(4*pi);

    
      end
  end
  Wr = 2* real(sum(segment_wr, 'all'));

end