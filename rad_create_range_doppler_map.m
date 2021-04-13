function result = rad_create_range_doppler_map(waveform,radar_platform,target_platform)
%RAD_CREATE_RANGE_DOPPLER_MAP Create a range-Doppler map
%   
%     result = NaN;
    
    radar_signal = create_signal(waveform);
    target_response = compute_target_response(radar_signal,radar_platform,target_platform);
    matched_filter_signal = compute_matched_filter_signal(radar_signal,target_response);
    doppler_processed_signal = compute_Doppler_processed_signal(radar_signal,matched_filter_signal);
    normalized_Doppler_signal = normalize_Doppler_processed_signal(doppler_processed_signal);
    noisy_range_Doppler_map = add_receiver_noise(normalized_Doppler_signal);
    
    result = noisy_range_Doppler_map;
    
end

function result = get_c()
    result = 3e8;
end

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

function result = compute_pulse_repetition_interval(pulse_repetition_frequency)
    result = 1/pulse_repetition_frequency;
end

function result = get_samples_during_pulse(num_fast_time_samples,...
    duty_factor)
    result = ceil(num_fast_time_samples*duty_factor);
end

function result = create_signal_envelope(prf,num_fast_time_samples,duty_factor)
    result.t = linspace(0,...
        compute_pulse_repetition_interval(prf),...
        num_fast_time_samples);
    result.voltage = zeros(1,num_fast_time_samples);
    
    samples_during_pulse = get_samples_during_pulse(num_fast_time_samples,...
        duty_factor);
    result.voltage(1:samples_during_pulse+1) = 1;
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
    target_slant_range = compute_target_slant_range(radar_platform.altitude,...
        target_platform.altitude,target_platform.range);
    target_doppler = compute_target_doppler(radar_signal.waveform.frequency,...
        radar_platform.velocity,target_platform.velocity);
    target_phase = compute_target_phase_for_each_pulse(radar_signal.waveform,...
        target_slant_range,target_doppler);
    time_shifted_signal = compute_time_shifted_signal(target_slant_range,...
        radar_signal);
    
    result = zeros(radar_signal.waveform.num_pulses,...
        radar_signal.waveform.num_fast_time_samples);
    result = result + target_phase*time_shifted_signal;
    result = blank_receiver_during_transmit(result,...
        get_samples_during_pulse(radar_signal.waveform.num_fast_time_samples,...
            radar_signal.waveform.duty_factor));
end

function result = compute_target_slant_range(radar_altitude,target_altitude,target_range)
    altitude_difference = radar_altitude - target_altitude;
    result = sqrt(altitude_difference^2 + target_range^2);
end

function result = compute_target_doppler(signal_frequency,radar_velocity,target_velocity)
    target_relative_velocity = compute_target_relative_velocity(radar_velocity,...
        target_velocity);
    radar_wavelength = compute_wavelength(signal_frequency);
    result = 2*target_relative_velocity/radar_wavelength;
end

function result = compute_target_relative_velocity(radar_velocity,target_velocity)
    result = radar_velocity - target_velocity;
end

function result = compute_target_phase_for_each_pulse(waveform,slant_range,target_doppler)
    pulse_times = compute_pulse_times(waveform);
    wave_number = compute_wave_number(waveform.frequency);
    result = exp(1j*2*pi*target_doppler*pulse_times)*exp(-1j*2*wave_number*slant_range);
end

function result = compute_time_shifted_signal(slant_range,radar_signal)
    unambiguous_range_swath = compute_unambigious_range(...
        radar_signal.waveform.pulse_repetition_frequency);
    time_delay = mod(slant_range,unambiguous_range_swath)/get_c();
    result = circshift(radar_signal.voltage,...
        ceil(time_delay*...
            radar_signal.waveform.num_fast_time_samples*...
            radar_signal.waveform.pulse_repetition_frequency));
end

function result = compute_unambigious_range(prf)
    result = 0.5*get_c()*compute_pulse_repetition_interval(prf);
end

function result = blank_receiver_during_transmit(reflected_signal,samples_during_pulse)
    result = reflected_signal;
    result(:,1:samples_during_pulse) = 0;
end

function result = compute_pulse_times(waveform)
    num_pulses = waveform.num_pulses;
    pulses = (0:num_pulses-1)';
    result = pulses*compute_pulse_repetition_interval(waveform.pulse_repetition_frequency);
end

function result = compute_wave_number(signal_frequency)
    radar_wavelength = compute_wavelength(signal_frequency);
    result = 2*pi/radar_wavelength;
end

function result = compute_wavelength(frequency)
    result = get_c()/frequency;
end

function result = compute_matched_filter_signal(radar_signal,target_response)
    result = zeros(radar_signal.waveform.num_pulses,...
        2*radar_signal.waveform.num_fast_time_samples-1);
    for i_pulse = 1:radar_signal.waveform.num_pulses
        result(i_pulse,:) = xcorr(target_response(i_pulse,:),...
            conj(radar_signal.voltage));
    end
    result(:,1:radar_signal.waveform.num_fast_time_samples-1) = [];
end

function result = compute_Doppler_processed_signal(radar_signal,matched_filter_signal)
    result = zeros(size(matched_filter_signal));
    for i_sample = 1:radar_signal.waveform.num_fast_time_samples
        result(:,i_sample) = fftshift(fft(matched_filter_signal(:,i_sample)));
    end
    
end

function result = normalize_Doppler_processed_signal(doppler_processed_signal)
    result = doppler_processed_signal./max(abs(doppler_processed_signal));

end

function result = add_receiver_noise(normalized_Doppler_processed_signal)
    peak_snr_dB = 10; % TODO: Refactor into input parameter or calculate from some radar system description
    peak_snr_linear = db2mag(peak_snr_dB);
    noise = sqrt(1/peak_snr_linear)*...
        (randn(size(normalized_Doppler_processed_signal))+...
         1j*randn(size(normalized_Doppler_processed_signal)));
    result = abs(normalized_Doppler_processed_signal + noise);
end

% function testcompute_pulse_repetition_interval(testCase)
%     assertEqual(testCase^-1,compute_pulse_repetition_interval(testCase));
% end

% 
% waveform = struct('frequency', 1.0e9,...
%                   'pulse_repetition_frequency', 50e3,...
%                   'duty_factor', 0.1,...
%                   'num_pulses', 64,...
%                   'num_fast_time_samples', 1000,...
%                   'bandwidth', 5e6,...
%                   'starting_phase', 0,...
%                   'chirp_direction', +1);
% radar_platform = struct('altitude',9144,...
%                         'velocity',300);
% target_platform = struct('range',10000,...
%                          'altitude',9144,...
%                          'velocity',-300);
