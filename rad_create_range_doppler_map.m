function result = rad_create_range_doppler_map(radar_signal,radar_platform,target_platform)
%RAD_CREATE_RANGE_DOPPLER_MAP Create a range-Doppler map
%   
%
    rng(0);
    target_response = compute_target_response(radar_signal,radar_platform,target_platform);
    matched_filter_signal = compute_matched_filter_signal(radar_signal,target_response);
    doppler_processed_signal = compute_doppler_processed_signal(radar_signal,matched_filter_signal);
    normalized_doppler_signal = normalize_doppler_processed_signal(doppler_processed_signal);
    noisy_range_doppler_map = add_receiver_noise(normalized_doppler_signal);
    
    result = noisy_range_doppler_map;
    
end


