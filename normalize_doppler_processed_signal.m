function result = normalize_doppler_processed_signal(doppler_processed_signal)
%NORMALIZE_DOPPLER_PROCESSED_SIGNAL
    result = doppler_processed_signal./max(abs(doppler_processed_signal(:)));

end

