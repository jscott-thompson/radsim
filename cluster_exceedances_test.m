rad_create_range_doppler_map_test;

vels = 0.5*compute_wavelength(waveform.frequency)*waveform.pulse_repetition_frequency*...
    linspace(-0.5,0.5,waveform.num_pulses);
detector_params = struct('num_training_cells',6,'num_guard_cells',3,...
    'cluster_distance_threshold',0.5,...
    'ranges',radar_signal.t*get_c(),...
    'vels',vels);
exceedances = detect_cfar_exceedances(range_doppler_map,detector_params);
detections = cluster_exceedances(exceedances,detector_params);

load gold_standard_detections.mat

assert(isequaln(detections,DETECTIONS));
