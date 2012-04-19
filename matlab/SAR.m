function R = SAR(s)

%% Parameters
global c delta_x delta_y f_carrier delta_f;
jj = 1j;
Dx = delta_x;
Dy = delta_y;
f0 = f_carrier;
Df = delta_f;
Z1 = -1;

%% Compute A
A = zeros(size(s));
for n = 1:size(s,3)
    A(:,:,n) = fftn(s(:,:,n));
end
Nx = size(A,1);
Ny = size(A,2);
N = size(A,3);
disp('FFT done');

%% Compute D
[i j n] = meshgrid(0:Nx-1, 0:Ny-1, 0:N-1);
  function Cijn = compute_C(i,j,n)
      
    if(i < Nx/2)
      kx = i/Nx * 2*pi/Dx;
    else
      kx = (i-Nx)/Nx * 2*pi/Dx;
    end
    
    if(j < Ny/2)
      ky = j/Ny * 2*pi/Dy;
    else
      ky = (j-Ny)/Ny * 2*pi/Dy;
    end

    k = 2*pi*(f0 + n*Df)/c;
    kz = sqrt(4*k^2 - kx^2 - ky^2);
    
    % Equation from Concealed Weapon paper
    Dijn = exp(-jj * kz * Z1);
    
    % Equation from SAR Seismic paper
    % t0 = abs(2*Z1/c);
    % Dijn = exp(-1j * (c*t0*k^2/2 - kz*abs(Z1))) * kz/k;
    % Dijn = exp(-1j * (- kz*abs(Z1)));% * abs(kz/k);
    
    Cijn = A(i+1,j+1,n+1) * Dijn;
  end
D = arrayfun(@compute_C, i, j, n);
disp('C computed');

%% Stolt Interpolation Indices
f_max = f0 + (N-1)*Df;
k_min = 2*pi*f0/c;
k_max = 2*pi*f_max/c;
kx_max = (Nx-1)/Nx * 2*pi/Dx / 2;
ky_max = (Ny-1)/Ny * 2*pi/Dy / 2;

kz_min = sqrt(4*k_min^2 - kx_max^2 - ky_max^2)
kz_max = sqrt(4*k_max^2 - 0 - 0)
kz = linspace(kz_min, kz_max, N);
    
%% Interpolation
n_valid_pts = zeros(Nx, Ny);
ns = zeros(Nx, N);
E = zeros(size(D));
for i = 0:Nx-1
  for j = 0:Ny-1
    
    % Interpolation Indices
    if(i < Nx/2)
      kx = i/Nx * 2*pi/Dx;
    else
      kx = (i-Nx)/Nx * 2*pi/Dx;
    end
    
    if(j < Ny/2)
      ky = j/Ny * 2*pi/Dy;
    else
      ky = (j-Ny)/Ny * 2*pi/Dy;
    end
    
    k = 0.5*sqrt(kx^2 + ky^2 + kz.^2);
    w = c*k;
    f = w/(2*pi);
    n = (c*k/(2*pi) - f0)/Df;

    data = D(i+1,j+1,:);
    E(i+1,j+1,:) = interp1( 1:N, data(:)', n + 1, 'linear', 0 );
  end
end

disp('E computed');

%% FFT
R = ifftn(E);
disp('R computed');

end