function calc_acf_and_msd(tr, n_tr, fps, delta_T, tau_max, min_t)    

tau_fr = 1:fps*tau_max; %tau in units of frames
tau = tau_fr/fps;%tau in units of seconds
dist = round(delta_T*fps);
min_dp = min_t*60*fps;
msd = zeros(n_tr, tau_fr(end));
msd_err = zeros(n_tr, tau_fr(end));
msd_std = zeros(n_tr, tau_fr(end));
acf = zeros(n_tr, tau_fr(end));
acf_err = zeros(n_tr, tau_fr(end));
acf_std = zeros(n_tr, tau_fr(end));
tr_len = zeros(n_tr, 1);

%set fit options for acf calculations    
my_acf = fittype('v_c^2*exp(-tau./tau_c)', 'independent', 'tau');
%intial values for non linear fit
tau_c0 = 1.5;
v_c0 = 1.5; 
opts_acf = fitoptions(my_acf);
set(opts_acf,'TolFun',1E-6, 'TolX', 1E-6, 'StartPoint', ...
    [tau_c0, v_c0]); 

%for msd calculations    
my_msd = fittype('2.*tau_c*(v_c^2).*(tau-tau_c.*(1-exp(-tau./tau_c)))+A',...
    'independent', 'tau');
A0 = 0.005;
tau_c0 = 1.5;
v_c0 = 1.5;
opts_msd = fitoptions(my_msd);
set(opts_msd,'TolFun',1E-6, 'TolX', 1E-6, 'StartPoint', ...
    [A0, tau_c0, v_c0]); 

n_ds = tr(end,1);
k=1;
for i=1:n_ds
    ind_ds = find(tr(:,1)==i);
    if isempty(ind_ds)
        continue;
    else
        tr_ds = tr(ind_ds, 2:end);
        tr_max = tr_ds(end,1);
        for j=1:tr_max
            ind_tr = find(tr_ds(:,1)==j);
            if length(ind_tr)>min_dp
                disp(['Track ',int2str(k),'/',int2str(n_tr),' processed!'])
                tr_len(k) = length(ind_tr)/fps;
                pos = tr_ds(ind_tr,3:4);
                time = tr_ds(ind_tr, 2)*60;
                %calculate correlation time and characteristic velocity via
                %msd 
                [msd(k,:), msd_err(k,:), msd_std(k,:)] = ...
                                                calc_msd(pos, tau_fr(end));
                %calc correlation time and characteristic velocity via 
                %angle correlation
                [acf(k,:), acf_err(k,:), acf_std(k,:)] = ...
                                    calc_acf(pos, time, max(tau_fr), dist);
                k = k + 1;
            end
        end
    end
end
msd = mean(msd)';
acf = mean(acf)';
msd_std = mean(msd_std)';
acf_std = mean(acf_std)';
msd_err = sqrt(mean(msd_err.^2))';
acf_err = sqrt(mean(acf_err.^2))';

%down sampling of data
% [msd, tmp] = down_sampling(msd, tau, fps, 0.5);
% [msd_err, tmp] = down_sampling_err(msd_err, tau, fps, 0.5);
[acf, tmp] = down_sampling(acf, tau, fps, 1);
[acf_err, tau_acf] = down_sampling_err(acf_err, tau, fps, 1);

acf_w = 1./(acf_err').^2;
msd_w = 1./(msd_err).^2;

assignin('base', 'tau', tau);
assignin('base', 'msd', msd);
assignin('base', 'tau_acf', tau_acf);
assignin('base', 'acf', acf);
assignin('base', 'msd_err', msd_err);
assignin('base', 'acf_err', acf_err);
assignin('base', 'msd_std', msd_std);
assignin('base', 'acf_std', acf_std);
assignin('base', 'msd_w', msd_w);
assignin('base', 'acf_w', acf_w);
assignin('base', 'tr_len', tr_len);

set(opts_acf, 'Weights', acf_w);
%do exponetial fit with my_acf;
if (sum(isnan(acf))==0)
    [acf_fit, acf_gof] = fit(tau_acf', acf', my_acf, opts_acf); 
    figure(1)
    dp = length(tau_acf);
    ind_ds = 1:1:dp;
    errorbar(tau_acf(ind_ds), acf(ind_ds), acf_err(ind_ds), '+r')
    hold on
    plot(acf_fit)
    hold off
    xlabel('Time \tau [s]')
    ylabel('ACF [µm^2/s^2]')
    xlim([tau_acf(1) tau_acf(end)]);
    ylim([min(acf) max(acf)]);
else
    acf_fit = NaN;
    acf_gof = NaN;
end
assignin('base', 'acf_fit', acf_fit);
assignin('base', 'acf_gof', acf_gof);

set(opts_msd, 'Weights', msd_w);

if (sum(isnan(msd))==0)
    [msd_fit, msd_gof] = fit(tau', msd,my_msd, opts_msd);
    figure(2)
    dp = length(tau);
    ind_ds = 1:1:dp;
    errorbar(tau(ind_ds), msd(ind_ds), msd_err(ind_ds), '+r')
    hold on
    plot(msd_fit)
    hold off
    xlabel('Time \tau [s]')
    ylabel('MSD [µm^2]')
    xlim([tau(1) tau(end)]);
    ylim([min(msd) max(msd)]);
else
    msd_fit = NaN;
    msd_gof = NaN;
end
assignin('base', 'msd_fit', msd_fit);
assignin('base', 'msd_gof', msd_gof);


function [acf, acf_err, acf_std] = calc_acf(pos, time, tau_max, dist)

tr_len = size(pos,1);
dp = tr_len - dist;
v = zeros(dp,2);
for j=1:dp;
    ind = j + dist;
    if (ind <= tr_len)
        s_x = pos(ind,1)-pos(j,1);
        s_y = pos(ind,2)-pos(j,2);
        delta_t = time(ind)-time(j);   
        v(j,1) = s_x/delta_t;
        v(j,2) = s_y/delta_t;
    end
end

acf = zeros(tau_max,1);
acf_err = zeros(tau_max,1);
acf_std = zeros(tau_max,1);

for j=1:tau_max;
    acf_j = zeros(dp-j,1);
    for l=1:dp-j
        ind = l+j;
        acf_j(l) = (v(l,1)*v(ind,1)+v(l,2)*v(ind,2));
    end
    if ~isempty(acf_j)
        acf(j) = mean(acf_j);
        acf_err(j) = std(acf_j)/sqrt(dp-j);
        acf_std(j) = std(acf_j);
    end
end

function [msd, msd_err, msd_std] = calc_msd(pos,tau_max)

tr_len = size(pos,1);
msd = zeros(tau_max,1);
msd_err = zeros(tau_max,1);
msd_std = zeros(tau_max,1);

for j=1:tau_max;
    sd = zeros(tr_len-j,1);
    for l=1:tr_len-j
        ind = l+j;
        s_vec = [pos(ind,1)-pos(l,1);pos(ind,2)-...
            pos(l,2)]; 
        sd(l) = norm(s_vec)^2;                                    
    end
    if ~isempty(sd)
        msd(j) = mean(sd);
        msd_err(j) = std(sd)/sqrt(tr_len-j);
        msd_std(j) = std(sd);
    end
end 

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

function [data_ds, time_ds] = down_sampling_err(data, time, f_s, f_ds)

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
    data_ds(i) = sqrt(mean(data((i-1)*delta_fr+1:i*delta_fr).^2));
    time_ds(i) = mean(time((i-1)*delta_fr+1:i*delta_fr));
end