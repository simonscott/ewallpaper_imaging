function [n_targets, targets_x, targets_y, targets_z] = generate_scene()
% Returns the x, y and z co-ordinates of the targets in the scene.
% The co-ordinates are in metres.
% The antenna array is assumed to lie in the XY plane, at z=0, centred
% around the z axis.

n_targets = 2;

targets_x = [0.5, -0.5];
targets_y = [0.5, -0.5];
targets_z = [2, 2];

end