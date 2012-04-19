function data = read_raw_file(filename)
Nx = 128;
Ny = 128;
Nz = 256;
data = zeros(Nx, Ny, Nz);
n = 1;

f = fopen(filename ,'r');
if f == -1
    disp('Cannot find specified file');
    return;
end

[values, count] = fscanf(f, '%f, %f\n');

for x = 1:Nx
  for y = 1:Ny
    for z = 1:Nz
      data(x, y, z) = values(n) + j*values(n+1);
      n = n + 2;
    end
  end
end

fclose(f);