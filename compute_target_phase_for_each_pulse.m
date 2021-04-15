function result = compute_target_phase_for_each_pulse(waveform,slant_range,target_doppler)
%COMPUTE_TARGET_PHASE_FOR_EACH_PULSE As target moves throughout the CPI
    pulse_times = compute_pulse_times(waveform);
    wave_number = compute_wave_number(waveform.frequency);
    result = exp(1j*2*pi*target_doppler*pulse_times)*...
             exp(-1j*2*wave_number*slant_range);
end

