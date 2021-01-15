function params = rad_load_run_matrix(run_matrix_file)
% RAD_LOAD_RUN_MATRIX Loads the run matrix from CSV file
%
% Usage: 
%   params = rad_load_run_matrix(run_matrix_file);
%
% Columns in the run matrix file should be defined in this order:
%  PRF
%  Duty_factor
%  Frequency
%  Radar_vel
%  Target_vel
%  Radar_alt
%  Target_alt
%  Target_rcs
%  Clutter_rcs_sd
%  Slant_range
%  N_pulses
%  N_targets
%  Target_vel2
%  Target_alt2
%  Target_rcs2
%  Slant_range2

  params = csvread(run_matrix_file,1,0);

  return

