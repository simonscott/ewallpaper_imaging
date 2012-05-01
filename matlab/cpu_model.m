function [total_time, comp_time, comm_time] = cpu_model(latency, bandwidth, ...
        f_add_cycles, f_multiply_cycles, sqrt_cycles, sin_cycles, clock_speed, mem_cycles, Nx, Ny, Nf,
      precompute_fft, precompute_phase_operator, precompute_interpolation)
  
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

  function time = fast_fft_time(N)
    % For a single FFT of length N:
    % [M] : (1/2)*N*log(N)*(4 LOADS and 2 STORES)
    % [Ops] : 8*N*log(N) real floating point operations
    mem_time = 0.5 * N * log2(N) * (do_mem(4) + do_mem(2));
    ops_time = N * log2(N) * (complex_add() + complex_multiply());
    time = mem_time + ops_time;
  end

  function time = slow_fft_time(N)
    coeff_time = (f_mult_cycles + 2 * sin_cycles) / clock_speed;
    mem_time = 0.5 * N * log2(N) * (do_mem(4) + do_mem(2));
    ops_time = N * log2(N) * (complex_add() + coeff_time * complex_multiply());
    time = mem_time + ops_time;    
  end

  function time = fft_time(N)
    if(precompute_fft)
      time = fast_fft_time(N);
    else
      time = slow_fft_time(N);
  end

  function time = fast_phase_operator()
    % Precompute phase operator phi (Nf complex values)
    % [M] : 2*Nf LOADS and Nf STORES
    % [Ops] : Nf Complex Multiplies
    mem_time = do_mem(2*Nf) + do_mem(Nf);
    ops_time = Nf * complex_multiply();
    time = mem_time + ops_time;
  end

  function time = slow_phase_operator()
    cycles = 3 * f_multiply_cycles + 2 * f_add_cycles + sqrt_cycles + 2 * sin_cycles;
    single_time = cycles / clock_speed + complex_multiply() + 2*do_mem(Nf);
    time = Nf * single_time;
  end
  
  function time = phase_operator()
    if(precompute_phase_operator)
      time = fast_phase_operator();
    else
      time = slow_phase_operator();
  end

  function time = fast_linear_interpolator()
    % [M] : Nf times : 2 COMPLEX LOADS, 1 REAL LOAD, 1 COMPLEX STORE
    % [Ops] Nf times :
    % 1 floor, 1 ceil
    % 2 real subtract
    % 2 complex.scalar multiplies, 1 complex add
    mem_time = Nf * do_mem(3.5);
    complex_scalar_time = 2 * f_multiply_cycles / clock_speed;
    single_ops_time = 4 * f_add_cycles / clock_speed + ...
               2 * complex_scalar_time + complex_add();
    time = Nf*single_ops_time + mem_time;
  end

  function time = slow_linear_interpolator()
    single_precompute_time = (3 * f_add_cycles + 2 * f_mult_cycles + sqrt_cycles) / clock_speed;
    mem_time = do_mem(3);
    complex_scalar_time = 2 * f_multiply_cycles / clock_speed;
    single_ops_time = 4 * f_add_cycles / clock_speed + ...
               2 * complex_scalar_time + complex_add();    
    time = Nf*(single_precompute_time + mem_time + single_tops_time);
  end
  
  function time = linear_interpolator()
    if(precompute_interpolator)
      time = fast_linear_interpolator();
    else
      time = slow_linear_interpolator();
  end

  % Initialize Timer
  total_time = 0;
  comp_time = 0;
  comm_time = 0;
  
  % Step 1: send_row
  comm_time = comm_time + send_line();
  
  % Step 2: (Nx - 1) receive and forward
  comm_time = comm_time + receive_and_forward(Nx-1, Nf);
  
  % Step 3: (2x) 1D FFT of length Nx
  comp_time = comp_time + 2*fft_time(Nx);
  
  % Step 4: send_col
  comm_time = comm_time + send_line();
  
  % Step 5: (Ny - 1) receive and forward
  comm_time = comm_time + receive_and_forward(Ny-1, Nf);
  
  % Step 6: (2x) 1D FFT of length Ny
  comp_time = comp_time + 2*fft_time(Ny);
  
  % Step 7: send_row
  comm_time = comm_time + send_line();
  
  % Step 8: (Nx - 1) receive and forward
  comm_time = comm_time + receive_and_forward(Nx-1, Nf);
  
  % Step 9: phase operator
  comp_time = comp_time + phase_operator();
  
  % Step 10: linear interpolator
  comp_time = comp_time + linear_interpolator();
  
  % Step 11: 1D IFFT of length Nf
  comp_time = comp_time + fft_time(Nf);
  
  % Step 12: send_col
  comm_time = comm_time + send_line();
  
  % Step 13: (Ny - 1) receive and forward
  comm_time = comm_time + receive_and_forward(Ny-1, Nf);
  
  % Step 14: (2x) IFFT of length Ny
  comp_time = comp_time + fft_time(Ny);
  
  % Step 15: send_row
  comm_time = comm_time + send_line();
  
  % Step 16: (Nx - 1) receive and forward
  comm_time = comm_time + receive_and_forward(Nx-1, Nf);
  
  % Step 17: (2x) IFFT of length Nx
  comp_time = comp_time + fft_time(Nx);
  
  % Step 18: send_col
  comm_time = comm_time + send_line();
  
  % Step 19: (Ny - 1) receive and forward
  comm_time = comm_time + receive_and_forward(Ny-1, Nf);
  
  total_time = comm_time + comp_time;
end