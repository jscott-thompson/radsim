frequency = 1e9;
radar_velocity = 300;
target_velocity = -1575;

result = compute_target_doppler(frequency,radar_velocity,target_velocity);


target_doppler = 2*(target_velocity-radar_velocity)*frequency/get_c();
assert(result == target_doppler);