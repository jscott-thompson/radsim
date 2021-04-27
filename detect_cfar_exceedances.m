function result = detect_cfar_exceedances(range_doppler_map,detector_params)
%DETECT_CFAR_EXCEEDANCES Detects cells that exceed CFAR threshold

[threshold,num_cells] = getCFARLevel(range_doppler_map,...
    detector_params.num_training_cells,...
    detector_params.num_guard_cells);
peak_snr = 10; peak_snr_linear = 10^(0.2*peak_snr);
sinr_ca_cfar = 10^(0.1*6); % Signal-to-interference-plus-noise ratio for cell-averaged CFAR required to achieve a certain Pfa/Pd spec
% prob_false_alarm = exp(-1*10*log10(sinr_ca_cfar)/peak_snr_linear);
detection_candidates = range_doppler_map./threshold>sinr_ca_cfar;
result = find(detection_candidates == 1);