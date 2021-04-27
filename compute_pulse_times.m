function result = compute_pulse_times(waveform)
%COMPUTE_PULSE_TIMES given constant PRF
    num_pulses = waveform.num_pulses;
    pulses = (0:num_pulses-1)';
    result = pulses*compute_pulse_repetition_interval(waveform.pulse_repetition_frequency);
end

