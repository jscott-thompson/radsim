function result = compute_time_shifted_signal(slant_range,radar_signal)
%COMPUTE_TIME_SHIFTED_SIGNAL to account for time delay of propagation
    unambiguous_range_swath = compute_unambiguous_range(...
        radar_signal.waveform.pulse_repetition_frequency);
    time_delay = mod(slant_range,unambiguous_range_swath)/get_c();
    result = circshift(radar_signal.voltage,...
        ceil(time_delay*...
            radar_signal.waveform.num_fast_time_samples*...
            radar_signal.waveform.pulse_repetition_frequency));
end

