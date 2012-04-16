% Declare all variable global
global c n_ant_x n_ant_y delta_x delta_y f_carrier bandwidth chirp_duration
global n_samps delta_f range_max

% Physical constants

c = 299792458;          % speed of light, in m/s

% Antenna array parameters

n_ant_x = 128;
n_ant_y = 128;
delta_x = 0.012;        % units are meters
delta_y = 0.012;        % units are meters

% RF system parameters

f_carrier = 10e9;       % units are Hz
bandwidth = 2e9;        % from 10 GHz to 12 GHz
chirp_duration = 4e-6;  % units are seconds
n_samps = 256;
delta_f = bandwidth / n_samps;

% Scene parameters

range_max = 6;          % maximum range to target, in metres
