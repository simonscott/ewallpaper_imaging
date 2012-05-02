% Script to call the cpu_model function

% Fixed parameters for Rocket processor
latency = 15e-9; % See Mullins (Design & Impl) and Hypertransport
f_add_cycles = 5;
f_mult_cycles = 5;
sqrt_cycles = 30; % ARM single precision square root
sin_cycles = 31; % ARM System Developer's Guide
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
                       sqrt_cycles, sin_cycles, ...
                       clock_speed, mem_cycles, Nx, Ny, Nf, ...
                       true, true, true);
                     
% comm_times = [];
% comp_times = [];
% % for bandwidth = 0.1e9 : 0.1e9 : 3e9
% for Nx = 16:256
%   Ny = Nx;
%   Nf = 2*Nx;
%   [total_time, comp_time, comm_time] = cpu_model(latency, bandwidth, ...
%                        f_add_cycles, f_mult_cycles, ... 
%                        sqrt_cycles, sin_cycles, ...
%                        clock_speed, mem_cycles, Nx, Ny, Nf, ...
%                        true, true, true);
%   comm_times(end+1,1) = comm_time;
%   comp_times(end+1,1) = comp_time;
% end

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