function result = rad_create_range_doppler_map(waveform,radar_platform,target_platform)
%RAD_CREATE_RANGE_DOPPLER_MAP Create a range-Doppler map
%   
%     result = NaN;
    
    radar_signal = create_signal(waveform);
    target_response = compute_target_response(radar_signal,radar_platform,target_platform);
    matched_filter_signal = compute_matched_filter_signal(target_response);
    doppler_processed_signal = compute_Doppler_processed_signal(matched_filter_signal);
    normalized_Doppler_signal = normalize_Doppler_processed_signal(doppler_processed_signal);
    noisy_range_Doppler_map = add_receiver_noise(normalized_Doppler_signal);
    
    result = noisy_range_Doppler_map;
    
end

function result = create_signal(waveform)
%CREATE_SIGNAL Create the radar signal to be transmitted
    signal_envelope = create_signal_envelope(waveform.pulse_repetition_frequency,...
        waveform.num_fast_time_samples_cpi,...
        waveform.duty_factor);
    result = add_lfm_baseband_modulation(signal_envelope,...
        compute_pulse_width(waveform.duty_factor,...
            compute_pulse_repetition_interval(waveform.pulse_repetition_frequency)),...
        waveform.bandwidth,...
        waveform.starting_phase,...
        waveform.chirp_direction);
end

function result = compute_pulse_repetition_interval(pulse_repetition_frequency)
    result = 1/pulse_repetition_frequency;
end

function result = compute_samples_single_pulse(num_fast_time_samples_cpi,...
    duty_factor)
    result = ceil(num_fast_time_samples_cpi*duty_factor);
end

function result = create_signal_envelope(prf,num_fast_time_samples_cpi,duty_factor)
    result.t = linspace(0,...
        compute_pulse_repetition_interval(prf),...
        num_fast_time_samples_cpi);
    result.voltage = zeros(1,num_fast_time_samples_cpi);
    
    samples_in_one_pulse = compute_samples_single_pulse(num_fast_time_samples_cpi,...
        duty_factor);
    result.voltage(1:samples_in_one_pulse+1) = 1;
end

function result = compute_pulse_width(duty_factor,pulse_repetition_interval)
    result = duty_factor*pulse_repetition_interval;
end

function result = add_lfm_baseband_modulation(signal_envelope,bandwidth,pulse_width,starting_phase,chirp_direction)
    result.t = signal_envelope.t;
    result.voltage = signal_envelope.voltage.*...
        exp(1j*chirp_direction*bandwidth/pulse_width*...
        (signal_envelope.t-starting_phase).^2);
end

function result = compute_target_response(radar_signal,radar_platform,target_platform)
    

end

function result = compute_matched_filter_signal(radar_signal,target_response)
    result = zeros(radar_signal.num_pulses_coherent_processing_interval,...
        2*radar_signal.num_fast_time_samples_cpi);
    for i_pulse = 1:radar_signal.num_pulses_coherent_processing_interval)
        result(i_pulse,:) = xcorr(target_response(i_pulse,:),...
            conj(radar_signal));
    end
    result(:,1:radar_signal.num_fast_time_samples_cpi-1) = [];
end

function result = compute_Doppler_processed_signal(matched_filter_signal)
    result = NaN;
    
end

function result = normalize_Doppler_processed_signal(doppler_processed_signal)
    result = NaN;

end

function result = add_receiver_noise(normalized_Doppler_processed_signal)
    result = NaN;

end

% function testcompute_pulse_repetition_interval(testCase)
%     assertEqual(testCase^-1,compute_pulse_repetition_interval(testCase));
% end

% 
% waveform = struct('frequency', 1.0e9,...
%                   'pulse_repetition_frequency', 50e3,...
%                   'duty_factor', 0.1,...
%                   'num_pulses_coherent_processing_interval', 64,...
%                   'num_fast_time_samples_cpi', 1000,...
%                   
%               
% basebandWaveform = struct('bandwidth', 5e6,...
%                           'phase_start', 0,...
%                           'chirp_direction', +1,...