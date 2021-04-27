test_compute_pulse_repetition_interval(1e9);

function test_compute_pulse_repetition_interval(testCase)
    assert(testCase^-1 == compute_pulse_repetition_interval(testCase));
end
