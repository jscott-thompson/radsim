% RAD_CREATE_LHS_RUN_MATRIX Creates a run matrix using LHS
%

input_range = [50000  100000;  % True target range (m)
                  50     300;  % True target radial velocity (m/s)
                  50     300;  % Radar platform velocity (m/s)
                1e-5    2e-5;  % PRI (s)
              2.5e-2 3.75e-2]; % Wavelength (m)
n_design_points = 100000;
n_design_iterations = 10;
X = lhsdesign(n_design_points,size(input_range,1),'criterion','maximin','Iterations',n_design_iterations);
% X = lhsdesign(1000,5);

scale_factor = (input_range(:,2)-input_range(:,1)) + input_range(:,1);
scaled_X = X.*scale_factor';

show_matrix_plot = true;
if show_matrix_plot
    plotmatrix(scaled_X);
end


show_iteration_plots = false;
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
output_run_matrix = [1./scaled_X(:,4), 0.5*ones(n,1), c./scaled_X(:,5),...
    scaled_X(:,3), scaled_X(:,2), 10000.*ones(n,1), 10000.*ones(n,1),...
    10.*ones(n,1), 50.*ones(n,1), scaled_X(:,1), 64.*ones(n,1),...
    zeros(n,1),zeros(n,1),zeros(n,1),zeros(n,1),zeros(n,1)];
csvwrite('stupidly_simple_run_matrix.csv',output_run_matrix);
    