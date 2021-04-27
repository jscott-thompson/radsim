function result = blank_receiver_during_transmit(reflected_signal,samples_during_pulse)
%BLANK_RECEIVER_DURING_TRANSMIT to remove signals that the radar can't hear
    result = reflected_signal;
    result(:,1:samples_during_pulse) = 0;
end

