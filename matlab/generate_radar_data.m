function s = generate_radar_data(scene_type)
% Simulates a radar imaging system, and returns the received ADC data.
% scene_type is either 'simple' or 'sphere'. Default is 'simple'.
% Returns a 3D matrix of the received samples.

if nargin == 0
    scene_type = 'simple';
end

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

s = zeros(n_ant_x, n_ant_y, n_samps);

% Generate scene
if strcmp(scene_type, 'simple')
    [n_targets, targets_x, targets_y, targets_z, targets_refl] = generate_scene_simple();
else
    [n_targets, targets_x, targets_y, targets_z, targets_refl] = generate_scene_sphere();
end    

% For each antenna in array
for ant_x_idx = 0 : n_ant_x - 1
        
    ant_x = (ant_x_idx - n_ant_x/2) * delta_x;
    fprintf('Antenna column %d of %d\n', ant_x_idx, n_ant_x);
    
    for ant_y_idx = 0 : n_ant_y - 1
    
        ant_y = (ant_y_idx - n_ant_y/2) * delta_y;
        
        % For each frequency step
        for f_idx = 0 : n_samps - 1
            
            f = f_carrier + f_idx * delta_f;
            rx = 0;
            
            % For each target
            for target = 1:n_targets
            
                % Calculate delay
                dist = sqrt( (targets_x(target) - ant_x)^2 + (targets_y(target) - ...
                            ant_y)^2 + (targets_z(target) - 0)^2 );
                time_delay = dist * 2 / c;
            
                phase_delay =  time_delay * f * 2*pi;
                
                % Add to received signal
                % rx = rx + 1 * cos(phase_delay) - 1j * cos(phase_delay - pi/2);
                rx = rx + targets_refl(target) * exp(-j * phase_delay);
            end

            % Add to response
            s(ant_x_idx+1, ant_y_idx+1, f_idx+1) = rx;
        end
    end
end

end


