function detected_targets = radsim(waveform,radar_platform,target_platform,detector_params)
%RADSIM Simulate a range-Doppler map and detect targets using cell-averaged
%Constant False Alarm Rate (CFAR) detection

radar_signal = create_signal(waveform);
range_doppler_map = rad_create_range_doppler_map(radar_signal,...
    radar_platform,target_platform);
detector_params.ranges = radar_signal.t*get_c();
detector_params.vels = 0.5*compute_wavelength(waveform.frequency)*waveform.pulse_repetition_frequency*...
    linspace(-0.5,0.5,waveform.num_pulses);
detected_targets = detect_targets(range_doppler_map,detector_params);