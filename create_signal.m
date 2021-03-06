function result = create_signal(waveform)
%CREATE_SIGNAL Create the radar signal to be transmitted
    signal_envelope = create_signal_envelope(waveform.pulse_repetition_frequency,...
        waveform.num_fast_time_samples,...
        waveform.duty_factor);
    result = add_lfm_baseband_modulation(signal_envelope,...
        compute_pulse_width(waveform.duty_factor,...
            compute_pulse_repetition_interval(waveform.pulse_repetition_frequency)),...
        waveform.bandwidth,...
        waveform.starting_phase,...
        waveform.chirp_direction);
    result.waveform = waveform;
end

