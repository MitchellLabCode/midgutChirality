function [total_twist, local_twist] = compute_twist(center_line, normal_collection)
    % Compute the twist of a rod using center_line and normal_collection.
    % Methods used in this function are adapted from method 2a Klenin & Langowski 2000 paper.

    % Parameters
    % ----------
    % center_line: 3D array (n_nodes, 3) containing float type data.
    %              Time history of rod node positions.
    % normal_collection: 3D array (n_elems, 3) containing float type data.
    %                    Time history of rod elements normal direction.
    %
    % Returns
    % -------
    % total_twist : float 
    % local_twist : N-2 x 1 float array
    %   

    [N, ~] = size(center_line);

    total_twist = 0;
    local_twist = zeros(N - 2, 1);  % Only consider interior nodes

    % Compute s vector

    % s is the tangent
    s = center_line(2:end, :) - center_line(1:end-1, :);
    
    % Compute tangents unit directors
    tangent = s ./ vecnorm(s, 2, 2);  % Equivalent to _batch_norm

    if size(tangent,1)==size(normal_collection,1)
        idx = 1:size(tangent,1);
    else
        idx = 1:(size(tangent,1)-1);
    end

    % Compute the projection of normal collection (d1) on normal-binormal plane
    projection_of_normal_collection = normal_collection(idx, :) - ...
        sum(tangent .* normal_collection(idx, :), 2) .* tangent;
    projection_of_normal_collection = projection_of_normal_collection ./ ...
        vecnorm(projection_of_normal_collection, 2, 2);

    % Eq 27 in Klenin & Langowski 2000
    % p is defined on interior nodes
    p = cross(s(1:end-1, :), s(2:end, :), 2);
    p = p ./ vecnorm(p, 2, 2);  % Equivalent to _batch_norm

    % Compute the angle we need to turn d1 around s to get p
    % Sign part tells whether d1 must be rotated ccw(+) or cw(-) around s
    alpha = sign(sum(cross(projection_of_normal_collection(1:end-1, :), p, 2) .* s(1:end-1, :), 2)) .* ...
        acos(sum(projection_of_normal_collection(1:end-1, :) .* p, 2));
    
    gamma = sign(sum(cross(p, projection_of_normal_collection(2:end, :), 2) .* s(2:end, :), 2)) .* ...
        acos(sum(projection_of_normal_collection(2:end, :) .* p, 2));

    % Normalize angles to a full rotation
    alpha = alpha / (2 * pi);
    gamma = gamma / (2 * pi);
    
    twist_temp = alpha + gamma;

    % Ensure twist is between (-1/2 to 1/2)
    twist_temp(twist_temp > 0.5) = twist_temp(twist_temp > 0.5) - 1;
    twist_temp(twist_temp < -0.5) = twist_temp(twist_temp < -0.5) + 1;

    % Check for NaN values (occur when tangents are parallel)
    twist_temp(isnan(twist_temp)) = 0;

    % Store twist values
    local_twist(:) = twist_temp;
    total_twist = sum(twist_temp);
end