function result = add_lfm_baseband_modulation(signal_envelope,bandwidth,pulse_width,starting_phase,chirp_direction)
%ADD_LFM_BASEBAND_MODULATION Superimpose a chirp over the signal
    result.t = signal_envelope.t;
    result.voltage = signal_envelope.voltage.*...
        exp(1j*chirp_direction*bandwidth/pulse_width*...
        (signal_envelope.t-starting_phase).^2);
end

