function [n_targets, targets_x, targets_y, targets_z] = generate_scene()
% Returns the x, y and z co-ordinates of the targets in the scene.
% The co-ordinates are in metres.
% The antenna array is assumed to lie in the XY plane, at z=0, centred
% around the z axis.

%n_targets = 3;

%targets_x = [0.5, -0.5, 0];
%targets_y = [0.5, -0.5, 0];
%targets_z = [1.5, 1.5, 1.8];

targets_x = [];
targets_y = [];
targets_z = [];
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
        n_targets = n_targets + 1;
        
        targets_x(1, end+1) = x;
        targets_y(1, end+1) = y;
        targets_z(1, end+1) = cz - z;
        n_targets = n_targets + 1;
    end
end
% center = [0 0 1];
% radius = 0.5;
% n_targets = 20;
% dx = 2*radius/sqrt(n_targets);
% dy = 2*radius/sqrt(n_targets);
% 
% targets_x = []; %zeros(1, n_targets);
% targets_y = []; %zeros(1, n_targets);
% targets_z = []; %zeros(1, n_targets);
% n_targets = 0;
% for x = center(1)-radius/2:dx:center(1)+radius/2
%     for y = center(2)-radius/2:dy:center(2)+radius/2
%         z = sqrt(radius^2 - (x - center(1))^2 - (y - center(2))^2) + center(3);
%         targets_x(1,end+1) = x;
%         targets_y(1,end+1) = y;
%         targets_z(1,end+1) = z;
%         n_targets = n_targets+1;
%     end
% end

%[x y] = meshgrid(-0.75:0.05:0.75, -0.75:0.05:0.75);
%z = x.^2 + y.^2

end