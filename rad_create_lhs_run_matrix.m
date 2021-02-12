% RAD_CREATE_LHS_RUN_MATRIX Creates a run matrix using LHS
%

input_range = [50000  100000;  % True target range (m)
                  50     300;  % True target radial velocity (m/s)
                  50     300;  % Radar platform velocity (m/s)
                1e-5    2e-5;  % PRI (s)
              2.5e-2 3.75e-2]; % Wavelength (m)

X = lhsdesign(1000,5,'criterion','maximin','Iterations',10000);
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
    