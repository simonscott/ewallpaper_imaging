function total_time = cpu_model(latency, bandwidth, f_add_cycles, f_multiply_cycles, ... 
                                clock_speed, mem_cycles, Nx, Ny, Nf)
  
  function time = do_mem(num_complex)
    time = num_complex * 2 * mem_cycles / clock_speed;
  end

  function time = send_line()
    % [M] : Nf STORES
    % Takes 4l time
    time = do_mem(Nf) + 4*latency;
  end

  function time = receive_and_forward(num_messages, Nf)
    % [M] : 2 LOADS, 2 STORES
    %   - swap send and receive pointers
    % A = message_size / bandwidth
    % l = latency
    % one cycle takes A + l time
    % thus total time = (Nx - 1) * (A + l)
    % ignoring memory op time.    
    A = Nf * 2 * 4 * 8 / bandwidth;
    time = num_messages * (A + latency);
  end

  function time = complex_add()
    time = 2 * f_add_cycles / clock_speed;
  end

  function time = complex_multiply()
    time = (4 * f_multiply_cycles + 2 * f_add_cycles) / clock_speed;
  end

  function time = fft_time(N)
    % For a single FFT of length N:
    % [M] : (1/2)*N*log(N)*(4 LOADS and 2 STORES)
    % [Ops] : 8*N*log(N) real floating point operations
    mem_time = 0.5 * N * log2(N) * (do_mem(4) + do_mem(2));
    ops_time = N * log2(N) * (complex_add() + complex_multiply());
    time = mem_time + ops_time;
  end

  function time = phase_operator()
    % Precompute phase operator phi (Nf complex values)
    % [M] : 2*Nf LOADS and Nf STORES
    % [Ops] : Nf Complex Multiplies
    mem_time = do_mem(2*Nf) + do_mem(Nf);
    ops_time = Nf * complex_multiply();
    time = mem_time + ops_time;
  end

  function time = linear_interpolator()
    % [M] : Nf times : 2 COMPLEX LOADS, 1 REAL LOAD, 1 COMPLEX STORE
    % [Ops] Nf times :
    % 1 floor, 1 ceil
    % 2 real subtract
    % 2 complex.scalar multiplies, 1 complex add
    mem_time = Nf * do_mem(3.5);
    complex_scalar_time = 2 * f_multiply_cycles / clock_speed;
    ops_time = 4 * f_add_cycles / clock_speed + ...
               2 * complex_scalar_time + complex_add();
    time = ops_time + mem_time;
  end

  % Initialize Timer
  total_time = 0;
  
  % Step 1: send_row
  total_time = total_time + send_line();
  
  % Step 2: (Nx - 1) receive and forward
  total_time = total_time + receive_and_forward(Nx-1, Nf);
  
  % Step 3: (2x) 1D FFT of length Nx
  total_time = total_time + 2*fft_time(Nx);
  
  % Step 4: send_col
  total_time = total_time + send_line();
  
  % Step 5: (Ny - 1) receive and forward
  total_time = total_time + receive_and_forward(Ny-1, Nf);
  
  % Step 6: (2x) 1D FFT of length Ny
  total_time = total_time + 2*fft_time(Ny);
  
  % Step 7: send_row
  total_time = total_time + send_line();
  
  % Step 8: (Nx - 1) receive and forward
  total_time = total_time + receive_and_forward(Nx-1, Nf);
  
  % Step 9: phase operator
  total_time = total_time + phase_operator();
  
  % Step 10: linear interpolator
  total_time = total_time + linear_interpolator();
  
  % Step 11: 1D IFFT of length Nf
  total_time = total_time + fft_time(Nf);
  
  % Step 12: send_col
  total_time = total_time + send_line();
  
  % Step 13: (Ny - 1) receive and forward
  total_time = total_time + receive_and_forward(Ny-1, Nf);
  
  % Step 14: (2x) IFFT of length Ny
  total_time = total_time + fft_time(Ny);
  
  % Step 15: send_row
  total_time = total_time + send_line();
  
  % Step 16: (Nx - 1) receive and forward
  total_time = total_time + receive_and_forward(Nx-1, Nf);
  
  % Step 17: (2x) IFFT of length Nx
  total_time = total_time + fft_time(Nx);
  
  % Step 18: send_col
  total_time = total_time + send_line();
  
  % Step 19: (Ny - 1) receive and forward
  total_time = total_time + receive_and_forward(Ny-1, Nf);
end