%% Test getCFARLevel
rng(0,'twister');
num_rows = 10;
num_cols = 10;
middle_row = floor(num_rows/2);
middle_col = floor(num_cols/2);

range_doppler_map = zeros(num_rows,num_cols) + rand(num_rows,num_cols);
range_doppler_map(middle_row,middle_col) = 10;

num_test_cells = 1;
num_half_guard_cells = 1;
[thresh,num_cells] = getCFARLevel(range_doppler_map,num_test_cells,num_half_guard_cells);
required_snr = 6;
detection_candidates = range_doppler_map > thresh*required_snr;


assert(find(detection_candidates) == 45);

load getCFARLevelTest_references.mat
assert(all(all(thresh == thresh_reference)));
assert(all(all(num_cells == num_cells_reference)));