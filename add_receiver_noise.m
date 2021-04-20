function result = add_receiver_noise(normalized_doppler_processed_signal)
%ADD_RECEIVER_NOISE Adds normally distributed noise, hard coded magnitude
    peak_snr_dB = 10; % TODO: Refactor into input parameter or calculate from some radar system description
    peak_snr_linear = 10^(0.2*peak_snr_dB); % TODO: Understand why the hard-coded 0.2, and refactor
    noise = sqrt(1/peak_snr_linear)*...
        (randn(size(normalized_doppler_processed_signal))+...
         1j*randn(size(normalized_doppler_processed_signal)));
    result = abs(normalized_doppler_processed_signal + noise);
end

