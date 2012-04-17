function display_3d(image)
% Makes a 3D plot of the provided image.
% 'image' is a 3D matrix
% A "dot" is plotted at each non-zero location in the 3D matrix.
% The colour of the dot is determined by the value at that location.

radar_params;

% [sp_x, sp_y, sp_z] = sphere();
% dim = size(image);
% figure();

% Make plot
% for x = 1:dim(1)
%     for y = 1:dim(2)
%         for z = 1:dim(3)
%             
%             % If non-zero, make a sphere at this point
%             if abs(image(x, y, z)) > 5
%                 colour = (abs(image(x, y, z)) * ones(size(sp_z)));
%                 surf_x = (sp_x/2 * abs(image(x, y, z)) / 10) + x;
%                 surf_y = (sp_y/2 * abs(image(x, y, z)) / 10) + y;
%                 surf_z = (sp_z/2 * abs(image(x, y, z)) / 10) + z;
%                 
%                 surf(surf_x, surf_y, surf_z, colour);
%                 hold('on');
%             end
%         end
%     end
% end

dim = size(image);
figure();

p = patch(isosurface(image, 2));
isonormals(image, p);
set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
alpha(0.5);
view(3);
camlight;
lighting gouraud

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


% Fix up aspect ratio and axes
axis([1 n_ant_x 1 n_ant_y 1 n_samps]);
daspect([1 1 1]);
end