function result = compute_target_relative_velocity(radar_velocity,target_velocity)
%COMPUTE_TARGET_RELATIVE_VELOCITY Given radar and target velocities
    result = target_velocity - radar_velocity;
end

