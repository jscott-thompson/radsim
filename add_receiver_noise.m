function result = add_receiver_noise(normalized_doppler_processed_signal)
%ADD_RECEIVER_NOISE Adds normally distributed noise, hard coded magnitude
    peak_snr_dB = 10; % TODO: Refactor into input parameter or calculate from some radar system description
    peak_snr_linear = db2mag(peak_snr_dB);
    noise = sqrt(1/peak_snr_linear)*...
        (randn(size(normalized_doppler_processed_signal))+...
         1j*randn(size(normalized_doppler_processed_signal)));
    result = abs(normalized_doppler_processed_signal + noise);
end

