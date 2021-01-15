% RAD_LOAD_CONSTANTS Loads constants for Rad
% 
%% Logical flags
plot_clutter_rd_map = true;
plot_target_rd_map = true;
shuffle_seed = false;
simulate_clutter_response = true;
simulate_receiver_noise = true;
simulate_target_response = true;

%%  constants and conversion factors 
r2d=180/pi;             %   radians to degrees 
d2r=pi/180;             %   degrees to radians 
Nmi2m=1852;             %   nautical miles to meters 
m2Nmi=1/1852;           %   meters to nautical miles 
f2m=0.3048;             %   feet to meters
c=2.998e8;              %   speed of light, m/sec 
Re=6300e3;              %   radius of earth, m

%% Lazy enum for input parameters
I_PRF = 1;
DUTY_FACTOR = 2;
FREQUENCY = 3;
RADAR_VEL = 4;
TARGET_VEL = 5;
RADAR_ALT = 6;
TARGET_ALT = 7;
TARGET_RCS = 8;
CLUTTER_RCS_SD = 9;
TARGET_RANGE = 10;
N_PULSES = 11;
N_TARGETS = 12;
TARGET_VEL2 = 13;
TARGET_ALT2 = 14;
TARGET_RCS2 = 15;
TARGET_RANGE2 = 16;

%% Lazy enum for results
N_TARGETS_DETECTED = 1;
RANGE = 2;
DOPPLER = 3;
RANGE2 = 4;
DOPPLER2 = 5;
RUN_TIME = 6;
