function [n_targets, targets_x, targets_y, targets_z] = generate_scene_sphere()
% Returns the x, y and z co-ordinates and reflectivity of the targets in the scene.
% The co-ordinates are in metres, and the target is a sphere.
% The antenna array is assumed to lie in the XY plane, at z=0, centred
% around the z axis.

targets_x = [];
targets_y = [];
targets_z = [];
targets_refl = [];
n_targets = 0;

cx = 0;
cy = 0;
cz = 3;
R_max = 0.5;

Rs = linspace(0, R_max, 12);
for Ri = 1:length(Rs)
    R = Rs(Ri);
    
    thetas = linspace(0, 2*pi, round(R/R_max * 16));
    for thetai = 1:length(thetas)
        theta = thetas(thetai);
        x = cx + R*cos(theta);
        y = cy + R*sin(theta);
        z = real(sqrt(R_max^2 - (x-cx)^2 - (y-cy)^2));
        
        targets_x(1, end+1) = x;
        targets_y(1, end+1) = y;
        targets_z(1, end+1) = cz + z;
        targets_refl(1, end+1) = 1;
        n_targets = n_targets + 1;
        
        targets_x(1, end+1) = x;
        targets_y(1, end+1) = y;
        targets_z(1, end+1) = cz - z;
        targets_refl(1, end+1) = 1;
        n_targets = n_targets + 1;
    end
end

end