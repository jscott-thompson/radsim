function result = compute_matched_filter_signal(radar_signal,target_response)
%COMPUTE_MATCHED_FILTER_SIGNAL Performs cross correlation on the received
%signal
    result = zeros(radar_signal.waveform.num_pulses,...
        2*radar_signal.waveform.num_fast_time_samples-1);
    for i_pulse = 1:radar_signal.waveform.num_pulses
        result(i_pulse,:) = xcorr(target_response(i_pulse,:),...
            conj(radar_signal.voltage));
    end
    result(:,1:radar_signal.waveform.num_fast_time_samples-1) = [];
end

