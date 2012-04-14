function display_3d(image)
% Makes a 3D plot of the provided image.
% 'image' is a 3D matrix
% A "dot" is plotted at each non-zero location in the 3D matrix.
% The colour of the dot is determined by the value at that location.

radar_params;

[sp_x, sp_y, sp_z] = sphere();
dim = size(image);
figure();

% Make plot
for x = 1:dim(1)
    for y = 1:dim(2)
        for z = 1:dim(3)
            
            % If non-zero, make a sphere at this point
            if abs(image(x, y, z)) ~= 0
                colour = abs(image(x, y, z) * ones(size(sp_z)));
                surf(sp_x/2 + x, sp_y/2 + y, sp_z/2 + z, colour);
                hold('on');
            end
        end
    end
end

% Fix up aspect ratio and axes
axis([0 n_ant_x 0 n_ant_y 0 n_samps 0 255]);
daspect(dim);

end