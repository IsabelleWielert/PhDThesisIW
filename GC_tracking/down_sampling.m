function [data_ds, time_ds] = down_sampling(data, time, f_s, f_ds)

if f_ds > f_s
    errordlg('Error in Downsampling! Downsampling frequency higher than sampling frequency')
    return;
end
delta_fr = round((f_s/f_ds));
num_dp_old = length(data);
num_dp = floor(num_dp_old/delta_fr);
data_ds = zeros(1, num_dp);
time_ds = zeros(1, num_dp);
for i=1:num_dp
    data_ds(i) = mean(data((i-1)*delta_fr+1:i*delta_fr));
    time_ds(i) = mean(time((i-1)*delta_fr+1:i*delta_fr));
end