% RAD_CREATE_LHS_RUN_MATRIX Creates a run matrix using LHS
%

input_range = [50000  100000;  % True target range (m)
                  50     300;  % True target radial velocity (m/s)
                  50     300;  % Radar platform velocity (m/s)
                1e-5    2e-5;  % PRI (s)
              2.5e-2 3.75e-2]; % Wavelength (m)

X = lhsdesign(1000,5);

scale_factor = (input_range(:,2)-input_range(:,1)) + input_range(:,1);
scaled_X = X.*scale_factor';
