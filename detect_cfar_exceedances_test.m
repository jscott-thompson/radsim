rad_create_range_doppler_map_test;

detector_params = struct('num_training_cells',6,'num_guard_cells',3);
exceedances = detect_cfar_exceedances(range_doppler_map,detector_params);

load gold_standard_cfar_detection_candidates.mat
assert(all(exceedances == DetCell));