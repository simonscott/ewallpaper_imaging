% Script to call the cpu_model function

% Fixed parameters for Rocket processor
latency = 15e-9; % See Mullins (Design & Impl) and Hypertransport
f_add_cycles = 5;
f_mult_cycles = 5;
clock_speed = 1e9;
mem_cycles = 2;

% Nominal values for variable parameters
bandwidth = 1e9;    % bits/s
Nx = 128;
Ny = 128;
Nf = 256;

% Compute model
[total_time, comp_time, comm_time] = cpu_model(latency, bandwidth, ...
                       f_add_cycles, f_mult_cycles, ... 
                       clock_speed, mem_cycles, Nx, Ny, Nf);

cpu_load = comp_time / total_time;
frame_rate = (1/total_time);

% Calculate computational time based on Hopper measurements
hopper_num_flops = 29252;
actual_comp_time = hopper_num_flops * (f_add_cycles + f_mult_cycles)/2 ...
                    / clock_speed;

% Display results
fprintf('Total time = %f, framerate = %f, CPU load = %f\n', ...
            total_time, frame_rate, cpu_load);
fprintf('Modeled computation time = %f, actual computation time = %f\n', ...
            comp_time, actual_comp_time);