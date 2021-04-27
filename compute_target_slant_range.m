function result = compute_target_slant_range(radar_altitude,target_altitude,target_range)
%COMPUTE_TARGET_SLANT_RANGE Given radar and target altitude and target
%ground range
    altitude_difference = radar_altitude - target_altitude;
    result = sqrt(altitude_difference^2 + target_range^2);
end