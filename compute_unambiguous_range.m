function result = compute_unambiguous_range(prf)
%COMPUTE_UNAMBIGUOUS_RANGE Given the pulse repetition frequency
    result = 0.5*get_c()*compute_pulse_repetition_interval(prf);
end

