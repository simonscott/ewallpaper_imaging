% Calculates the total communication time if data is sampled at all 128*128
% antennas, and then sent to a subset of those antennas for computation.
% Once the data is in the subset array, the communication time is the same
% as for the original array. This is because the number of hops decreases
% by X, but the packet sizes increase by X.
% However, the cost of moving the data to the subset cannot be ignored, and
% this is the reason why this whole scheme is just a bad idea.
%
% Plot shows width of processing array on X axis, and total comm time on y
% axis.


pkt_size = 256 * 2 * 4;
link_speed = 1e9;
bw_delay = (pkt_size * 8) / link_speed;
latency = 4e-9;
NUM_CYCLES = 6;

array_size = [1:128];

init_time = (bw_delay + latency) .* (127-array_size) + (bw_delay .* (128./array_size) + latency) .* (127-array_size);
comp_time = (NUM_CYCLES) * (bw_delay * (128./array_size) + latency) .* array_size;
time = init_time + comp_time;

plot(time)