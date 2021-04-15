function result = compute_wave_number(signal_frequency)
%COMPUTE_WAVE_NUMBER given the signal frequency
    radar_wavelength = compute_wavelength(signal_frequency);
    result = 2*pi/radar_wavelength;
end

