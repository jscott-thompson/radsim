function radsim_ml_datagen(run_matrix_file,results_file)
%RADSIM_ML_DATAGEN Generate a dataset for training statistical machine
%learning algorithms

params = rad_load_run_matrix(run_matrix_file);

N_runs = size(params,1);

show_waitbar = true;
if show_waitbar
    f_waitbar = waitbar(0,['Case 0 of ' num2str(N_runs)]);
end

runs = 1:N_runs;
run_time = zeros(N_runs,1);
total_run_timer = tic;

for i_run = runs
    if show_waitbar
        h_waitbar = waitbar(i_run/N_runs,f_waitbar,['Case ' num2str(i_run) ' of ' num2str(N_runs) ': Last run duration was ' num2str(run_time) ' s; ETC: ' num2str(N_runs*toc(total_run_timer)/i_run/60) ' min']);
    end
    this_run_timer = tic;


    waveform = struct('frequency', params(i_run,3),...
                      'pulse_repetition_frequency',params(i_run,1),...
                      'duty_factor', params(i_run,2),...
                      'num_pulses', params(i_run,11),...
                      'num_fast_time_samples', 1000,...
                      'bandwidth',5e6,...  % TODO: Make this a function of carrier frequency
                      'starting_phase', 0,...
                      'chirp_direction', +1);
     radar_platform = struct('altitude', params(i_run,6),...
                             'velocity', params(i_run,4));
     target_platform = struct('range',params(i_run,10),...
                              'altitude', params(i_run,7),...
                              'velocity', params(i_run,5),...
                              'rcs',params(i_run,8));
     detector_params = struct('num_training_cells',6,...
                              'num_guard_cells',3,...
                              'cluster_distance_threshold',0.5);

     detected_targets(i_run) = radsim(waveform,radar_platform,...
         target_platform,detector_params);
     
     run_time(i_run) = toc(this_run_timer);
end
save_detected_targets_to_csv(detected_targets,results_file);
total_run_time = toc(total_run_timer);

% TODO: write a function to save detection  results to a data file
% TODO: constrain detect_targets to detect exactly 1 target every time.
% What is a heuristic for a detection when CFAR doesn't find any?
