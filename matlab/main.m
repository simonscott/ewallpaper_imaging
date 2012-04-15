%s = generate_radar_data();
%save('scene_1.mat', 'scene_1');
load('scene_1.mat', 'scene_1');
image = SAR(scene_1);
display_3d(image);