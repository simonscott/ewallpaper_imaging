function R = SAR_2D(s)

%% Parameters
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

% Rename parameters

jj = 1j;
Dx = delta_x;
Dy = delta_y;
f0 = f_carrier;
Df = delta_f;
Z1 = -1;

A = fftn(s);
Nx = size(A,1);
Ny = size(A,2);
for i = 1:Nx-1
  for j = 1:Ny-1
    k = 2*pi*f0/c;
    kx = i/Nx * 2*pi/Dx;
    ky = j/Ny * 2*pi/Dy;
    kz = sqrt(4*k^2 - kx^2 - ky^2);
    A = A * exp(-j * kz * Z1);
  end
end

R = ifftn(A);