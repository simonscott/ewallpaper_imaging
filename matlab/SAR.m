function SAR

%% Parameters
jj = j;
Dx;
Dy;
f0;
Df;
c;
Z1;

%% Compute A
A = fftn(s);
Nx = size(A,1);
Ny = size(A,2);
N = size(A,3);

%% Compute D
[i j n] = meshgrid(1:Nx, 1:Ny, 1:N);
  function Cijn = compute_C(i,j,n)
    kz = 4*(2*pi*(f0 + n*Df)/c)^2 - (i/Nx * 2*pi/Dx)^2 - (j/Ny * 2*pi/Dy)^2;
    kz = sqrt(kz);
    Dijn = exp(-jj * kz * Z1);
    Cijn = A(i,j,n) * Dijn;
  end
D = arrayfun(@Cijn, i, j, n);

%% Stolt Interpolation Indices
kz_min = 2/c * 2*pi*f0;
kz_max = sqrt(4*(2*pi*(f0 + (N-1)*Df)/c)^2 - ((Nx-1)/Nx*2*pi/Dx)^2 - ((Ny-1)/Ny*2*pi/Dy)^2);
kz = linspace(kz_min, kz_max, N);

%% Interpolation
E = zeros(size(D));
  function interpolate(i,j)
    % Interpolation Indices
    n = 1/Df * (c/(2*pi)*sqrt((kz.^2 + (i/Nx*2*pi/Dx)^2 + (j/Ny*2*pi/Dy)^2)/4) - f0);
    E(i,j,:) = interp1(D(i,j,:), n);
  end
arrayfun(@interpolate, i, j);

%% FFT
R = ifftn(E);

end