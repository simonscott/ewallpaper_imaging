function R = SAR_simon(s)

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
[i j n] = meshgrid(1:Nx, 1:Ny, 1:N);
  function Cijn = compute_C(i,j,n)
      
    if(i <= Nx/2)
      kx = i/Nx * 2*pi/Dx;
    else
      kx = (i-Nx)/Nx * 2*pi/Dx;
    end
    
    if(j <= Ny/2)
      ky = j/Ny * 2*pi/Dy;
    else
      ky = (j-Ny)/Ny * 2*pi/Dy;
    end

    k = 2*(2*pi*(f0 + n*Df)/c);
    kz = sqrt(k^2 - kx^2 - ky^2);
    Dijn = exp(-jj * kz * Z1);
    Cijn = A(i,j,n) * Dijn;
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
for i = 1:Nx
  for j = 1:Ny
    
    fprintf('i=%d, j=%d\n', i ,j);
      
    % Interpolation Indices
    if(i <= Nx/2)
      kx = i/Nx * 2*pi/Dx;
    else
      kx = (i-Nx)/Nx * 2*pi/Dx;
    end
    
    if(j <= Ny/2)
      ky = j/Ny * 2*pi/Dy;
    else
      ky = (j-Ny)/Ny * 2*pi/Dy;
    end
    
    k = 0.5*sqrt(kx^2 + ky^2 + kz.^2);
    w = c*k;
    f = w/(2*pi);
    n = (c*k/(2*pi) - f0)/Df + 1;

    data = D(i,j,:);
    E(i,j,:) = interp1( 1:N, data(:)', n, 'linear', 0 );
  end
end

disp('E computed');

%% FFT
R = ifftn(E);
disp('R computed');

end