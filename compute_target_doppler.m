function result = compute_target_doppler(signal_frequency,radar_velocity,target_velocity)
%COMPUTE_TARGET_DOPPLER Given the signal freq, radar vel, and target vel
    target_relative_velocity = compute_target_relative_velocity(radar_velocity,...
        target_velocity);
    radar_wavelength = compute_wavelength(signal_frequency);
    result = 2*target_relative_velocity/radar_wavelength;
end

