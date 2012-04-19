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
data = values(1:2:end) + j*values(2:2:end);
data = reshape(data, [Nz Ny Nx]);
data = permute(data, [3 2 1]);

fclose(f);