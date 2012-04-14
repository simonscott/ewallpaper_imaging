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

% Rename parameters

jj = 1j;
Dx = delta_x;
Dy = delta_y;
f0 = f_carrier;
z0 = -2;

% Fake data
s_new = zeros(n_ant_x, n_ant_y);
for x = 1:n_ant_x
  x_pos = (x-n_ant_x/2) * delta_x;
  
  for y = 1:n_ant_y
    y_pos = (y-n_ant_y/2) * delta_y;
    dist = sqrt(x_pos^2 + y_pos^2 + z0^2);
    delta_t = 2*dist/c;
    delta_phase = delta_t * 2*pi*f0;

    dist2 = sqrt((x_pos-0.2)^2 + y_pos^2 + z0^2);
    delta_t2 = 2*dist2/c;
    delta_phase2 = delta_t2 * 2*pi*f0;
    
    s_new(x,y) = exp(-jj * delta_phase) + exp(-jj * delta_phase2);
  end
end

A = fftn(s);
Nx = size(A,1);
Ny = size(A,2);

for i = 0:Nx-1
  for j = 0:Ny-1
    k = 2*pi*f0/c;
    
    %kx = (i - Nx/2)/Nx * 2*pi/Dx;
    if(i <= Nx/2)
      kx = i/Nx * 2*pi/Dx;
    else
      kx = (i-Nx)/Nx * 2*pi/Dx;
    end
    
    %ky = (j - Ny/2)/Ny * 2*pi/Dy;
    if(j <= Ny/2)
      ky = j/Ny * 2*pi/Dy;
    else
      ky = (j-Ny)/Ny * 2*pi/Dy;
    end
    
    kz = sqrt(4*k^2 - kx^2 - ky^2);
    A(i+1,j+1) = A(i+1,j+1) * exp(-jj * kz * z0);
  end
end

R = ifftn(A);