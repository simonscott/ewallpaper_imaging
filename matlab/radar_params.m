% Physical constants

c = 299792458;      % speed of light, in m/s

% Antenna array parameters

n_ant_x = 128;
n_ant_y = 128;
delta_x = 0.025;    % units are meters
delta_y = 0.025;    % units are meters

% RF system parameters

f_carrier = 20e9;   % units are Hz
bandwidth = 1e9;    % units are Hz
n_samps = 128;
delta_f = bandwidth / n_samps;