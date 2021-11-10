function chi_d = guidance(position, waypoints, threshold_radius, lookahead_delta)
 
    % Find the two waypoints we are between
    distances_to_waypoints = vecnorm(waypoints - position);
    [~, closest_waypoint_indices] = sort(distances_to_waypoints);
    closest_waypoint_indices = sort(closest_waypoint_indices(1:2));
 
    % Check if we are close enough to the second closest way point to
    % swap to the next pair
    if distances_to_waypoints(closest_waypoint_indices(2)) < threshold_radius
        closest_waypoint_indices = closest_waypoint_indices + 1;
    end
 
    closest_waypoints = waypoints(:, closest_waypoint_indices);
 
    % Find the angle between the path vector and north axis
    path_vector = closest_waypoints(:, 2) - closest_waypoints(:, 1);
    pi_p = atan2(path_vector(2), path_vector(1)); 
 
    % Find the cross track error
    y_e_p = crosstrackWpt(closest_waypoints(1, 2),...
                          closest_waypoints(2, 2),...
                          closest_waypoints(1, 1),...
                          closest_waypoints(2, 1),...
                          position(1),...
                          position(2));
 
 
    Kp = 1/lookahead_delta;
    chi_d = pi_p - atan(Kp*y_e_p);
 
end