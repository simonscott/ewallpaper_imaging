function write_raw_file(data, filename)
f = fopen(filename ,'w');
Nx = 128;
Ny = 128;
Nz = 256;
for i = 1:Nx
  for j = 1:Ny
    for k = 1:Nz
      fprintf(f, '%f, %f\n', real(data(i,j,k)), imag(data(i,j,k)));
    end
  end
end
fclose(f);