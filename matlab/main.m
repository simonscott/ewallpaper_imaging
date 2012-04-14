%s = generate_radar_data();
%save('radar_raw_data.mat', 's');
load('radar_raw_data.mat', 's');
image = SAR(s);
display_3d(image);