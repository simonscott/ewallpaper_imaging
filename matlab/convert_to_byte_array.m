function data = convert_to_byte_array(data, filename)
   data = abs(data);
   minimum = min(data(:));
   maximum = max(data(:));
   data = (data - minimum)*255/(maximum - minimum);
   file = fopen(filename, 'w');
   fwrite(file, data, 'uint8');
   fclose(file);
end