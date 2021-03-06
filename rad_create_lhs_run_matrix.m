% RAD_CREATE_LHS_RUN_MATRIX Creates a run matrix using LHS
%

% input_range = [ 5000  100000;  % True target range (m)
%                 -300     300;  % True target radial velocity (m/s)
%                   50     300;  % Radar platform velocity (m/s)
%                 1e-5    2e-5;  % PRI (s)
%               2.5e-2 3.75e-2]; % Wavelength (m)
% n_design_points = 10000;
input_range = [ 5000  100000; % True target range (m)
                 300     300; % True taget radial velocity (m/s)
                  50      50; % Radar platform velocity (m/s)
                1e-3    1e-3; % PRI (s)
              0.0909  0.0909];% Wavelength (m)
n_design_points = 100;
n_design_iterations = 10;
X = lhsdesign(n_design_points,size(input_range,1),'criterion','maximin','Iterations',n_design_iterations);

scale_factor = (input_range(:,2)-input_range(:,1));
scaled_X = X.*scale_factor' + input_range(:,1)';

show_matrix_plot = false;
if show_matrix_plot
    plotmatrix(scaled_X);
end


show_iteration_plots = true;
if show_iteration_plots
    subplot(2,3,1);
    plot(scaled_X(:,1),'.');
    ylabel('Target Range (m)');
    
    subplot(2,3,2);
    plot(scaled_X(:,2),'.');
    ylabel('Target Rdot (m/s)');
    
    subplot(2,3,3);
    plot(scaled_X(:,3),'.');
    ylabel('Radar velocity (m/s)');
    
    subplot(2,3,4);
    plot(scaled_X(:,4),'.');
    ylabel('PRI (s)');
    
    subplot(2,3,5);
    plot(scaled_X(:,5),'.');
    ylabel('Wavelength (m)');
end

% TODO: Translate scaled design matrix into expected format for
% simpleRadarModel3
% Columns in the run matrix file should be defined in this order:
%  PRF                      1/PRI (scaled_X(:,4))
%  Duty_factor              Fixed 
%  Frequency                c/wavelength (scaled_X(:,5))
%  Radar_vel                (scaled_X(:,3))
%  Target_vel               (scaled_X(:,2))
%  Radar_alt                Fixed
%  Target_alt               Fixed
%  Target_rcs               Fixed
%  Clutter_rcs_sd           Fixed
%  Slant_range              (scaled_X(:,1))
%  N_pulses                 Fixed
%  N_targets                Fixed
%  Target_vel2              N/A
%  Target_alt2              N/A
%  Target_rcs2              N/A
%  Slant_range2             N/A
n = size(scaled_X(:,1),1);
i_duty_factor       = 0.5;
i_radar_alt         = 10000;
i_target_alt        = 10000;
i_target_rcs        = 10;
i_clutter_rcs_sd    = 200;
i_n_pulses          = 64;
i_n_targets         = 1;
% i_target_vel2              N/A
% i_target_alt2              N/A
% i_target_rcs2              N/A
% i_slant_range2             N/A

output_run_matrix = [1./scaled_X(:,4), i_duty_factor*ones(n,1), c./scaled_X(:,5),...
    scaled_X(:,3), scaled_X(:,2), i_radar_alt.*ones(n,1), i_target_alt.*ones(n,1),...
    i_target_rcs.*ones(n,1), i_clutter_rcs_sd.*ones(n,1), scaled_X(:,1), i_n_pulses.*ones(n,1),...
    i_n_targets*ones(n,1),zeros(n,1),zeros(n,1),zeros(n,1),zeros(n,1)];
csvwrite('stupidly_simple_run_matrix.csv',output_run_matrix);
    