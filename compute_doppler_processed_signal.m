function result = compute_doppler_processed_signal(radar_signal,matched_filter_signal)
%COMPUTE_DOPPLER_PROCESSED_SIGNAL by performing an FFT on the match
%filtered signal
    result = zeros(size(matched_filter_signal));
    for i_sample = 1:radar_signal.waveform.num_fast_time_samples
        result(:,i_sample) = fftshift(fft(matched_filter_signal(:,i_sample)));
    end
    
end

