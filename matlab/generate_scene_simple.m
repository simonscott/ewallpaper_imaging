function [n_targets, targets_x, targets_y, targets_z, targets_refl] = generate_scene_simple()
% Returns the x, y and z co-ordinates and the reflectivity of the targets in the scene.
% The co-ordinates are in metres.
% The antenna array is assumed to lie in the XY plane, at z=0, centred
% around the z axis.

n_targets = 3;

targets_x = [0.5, -0.5, 0];
targets_y = [0.5, -0.5, 0];
targets_z = [1.5, 1.5, 1.8];

targets_refl = [1, 1, 2];

end