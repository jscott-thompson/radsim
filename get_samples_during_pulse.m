function result = get_samples_during_pulse(num_fast_time_samples,...
    duty_factor)
%GET_SAMPLES_DURING_PULSE Gets the fast time samples during pulse transmit
%
    result = ceil(num_fast_time_samples*duty_factor);
end

