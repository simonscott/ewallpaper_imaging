function s = generate_radar_data()
%Simulates a radar imaging system, and returns the received ADC data.
%   Returns a 3D matrix of the received samples.

% Import system parameters
radar_params

% Ensure that Nyquist sampling requirements are met
max_delta_x = (c / (f_carrier + bandwidth)) / 2;
min_n_samps = (2 * range_max) / (c / (2*bandwidth));

if(delta_x >= max_delta_x)
    disp('Error: delta_x too large.');
end

if(n_samps <= min_n_samps)
    disp('Error: n_samps too small.');
end

% Calculate acheivable resolution
angle_subtended = (atan(n_ant_x * delta_x / range_max));
res_cross_range = (c/(f_carrier + bandwidth/2)) / (4*sin(angle_subtended/2));
res_range = c / (2*bandwidth);

out_str = sprintf('Range resolution = %0.2f m. Cross-range resolution = %0.2f m', ...
                res_range, res_cross_range);
disp(out_str);

% Discrete time for simulation of analogue circuitry
% The period of a 12 GHz wave is 83ps
timestep = 5e-12;

% Transmitted chirp signal
tx_time = 0 : timestep : chirp_duration;
tx = chirp(tx_time, f_carrier, chirp_duration, f_carrier + bandwidth);

% Plot FFT of chirp
fft_x_axis = 0 : 1 / chirp_duration : 1/timestep;
plot(fft_x_axis, abs(fft(tx)));

% Generate scene
[n_targets, targets_x, targets_y, targets_z] = generate_scene();

% For each antenna in array
for ant_x = 1:n_ant_x
    for ant_y = 1:n_ant_y
        
        % Send chirp
        
        % Record response 
        
        % Take FFT of response
        
    end
end

s = 0;

end


