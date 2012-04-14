function R = SAR(s)

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

% Scene parameters

range_max = 6;          % maximum range to target, in metres

jj = 1j;
Dx = delta_x;
Dy = delta_y;
f0 = f_carrier;
Df = delta_f;
Z1 = -1;

%% Compute A
%A = fftn(s);
A = zeros(size(s));
for n = 1:size(s,3)
    A(:,:,n) = fftn(s(:,:,n));
end
Nx = size(A,1);
Ny = size(A,2);
N = size(A,3);
disp('FFT done');

%% Compute D
[i j n] = meshgrid(1:Nx, 1:Ny, 1:N);
  function Cijn = compute_C(i,j,n)
    kz = 4*(2*pi*(f0 + n*Df)/c)^2 - (i/Nx * 2*pi/Dx)^2 - (j/Ny * 2*pi/Dy)^2;
    kz = sqrt(kz);
    Dijn = exp(-jj * kz * Z1);
    Cijn = A(i,j,n) * Dijn;
  end
D = arrayfun(@compute_C, i, j, n);
disp('C computed');

% %% Stolt Interpolation Indices
% kz_min = 2/c * 2*pi*f0;
% f_max = f0 + (N-1)*Df;
% k_max = 2*pi*f_max/c
% kx_max = (Nx-1)/Nx * 2*pi/Dx
% ky_max = (Ny-1)/Ny * 2*pi/Dy
% 
% kz_max = sqrt((2*k_max)^2 - kx_max^2 - ky_max^2)
% kz = linspace(kz_min, kz_max, N);
% kz
% 
% %% Interpolation
% E = zeros(size(D));
%   function interpolate(i,j)
%     % Interpolation Indices
%     term_1 = i/Nx * 2*pi/Dx;
%     term_2 = j/Ny * 2*pi/Dy;
%     k = 0.5*sqrt(term_1^2 + term_2^2 + kz.^2)
%     
%     w = c*k;
%     f = w/(2*pi)
%     
%     n = (c*k/(2*pi) - f0)/Df
%     
%     % n = 1/Df * (c/(2*pi)*sqrt((kz.^2 + (i/Nx*2*pi/Dx)^2 + (j/Ny*2*pi/Dy)^2)/4) - f0);
%     E(i,j,:) = interp1(D(i,j,:), n);
%   end
% arrayfun(@interpolate, i, j);
% disp('E computed');

%% FFT
R = ifftn(D);
disp('R computed');

end