function result = create_signal_envelope(prf,num_fast_time_samples,duty_factor)
%CREATE_SIGNAL_ENVELOPE Creates a binary array with 1s during transmit
    result.t = linspace(0,...
        compute_pulse_repetition_interval(prf),...
        num_fast_time_samples);
    result.voltage = zeros(1,num_fast_time_samples);
    
    samples_during_pulse = get_samples_during_pulse(num_fast_time_samples,...
        duty_factor);
    result.voltage(1:samples_during_pulse+1) = 1;
end
