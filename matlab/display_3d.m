function display_3d(image, display_method, threshold)
% Makes a 3D plot of the provided image.
% 'image' is a 3D matrix
% 'display_method' (optional) is either 'isosurf' (default) or 'sphere'.
% 'threshold' (optional) is the cutoff threshold for plotting.
% A "dot" is plotted at each non-zero location in the 3D matrix.
% The colour of the dot is determined by the value at that location.

radar_params;

if nargin == 1
    display_method = 'isosurf';
    threshold = 2;
end

if nargin == 2
    threshold = 2;
end

dim = size(image);
figure();

% Display data as a bunch of spheres at all non-zero points
% Size and color of sphere scales with signal magnitude.
if strcmp(display_method, 'sphere')
    
    [sp_x, sp_y, sp_z] = sphere();

    for x = 1:dim(1)
        for y = 1:dim(2)
            for z = 1:dim(3)

                % If non-zero, make a sphere at this point
                if abs(image(x, y, z)) > threshold
                    colour = (abs(image(x, y, z)) * ones(size(sp_z)));
                    surf_x = (sp_x/2 * abs(image(x, y, z)) / 10) + x;
                    surf_y = (sp_y/2 * abs(image(x, y, z)) / 10) + y;
                    surf_z = (sp_z/2 * abs(image(x, y, z)) / 10) + z;

                    surf(surf_x, surf_y, surf_z, colour);
                    hold('on');
                end
            end
        end
    end

% Display data as a bunch of surfaces connecting all points that have
% magnitude equal to the cutoff threshold
elseif strcmp(display_method, 'isosurf')

    p = patch(isosurface(image, threshold));
    isonormals(image, p);
    set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
    alpha(0.5);
    view(3);
    camlight;
    lighting gouraud
end

% Fix up aspect ratio and axes
axis([1 n_ant_x 1 n_ant_y 1 n_samps]);
daspect([1 1 1]);

% This doesn't work
% ss=['y','m','c','r','g','b','w','k','r','r'];
% k=1;
% for i=28:-3:1
%  p=patch(isosurface(image,i));
%  isonormals(image,p);
%  alpha(0.8);
%  hold on;
% 
%  set(p,'FaceColor',ss(k),'EdgeColor','none');
%  view(3);
%  camlight 
%  lighting gouraud
%  k=k+1;
% end

end