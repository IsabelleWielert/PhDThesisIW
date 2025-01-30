function eval_data


%global switching window in units of minutes
min_t = 19.5;
max_t = 21;


%open dialog for choosing txt.file
[FileName,PathName] = uigetfile('*.txt',...
                                'Choose text file containing speed data!');
%if cancel is pushed, then return.
if (isscalar(FileName) == 1) && (isscalar(PathName) == 1);
    return;
end  
cd(PathName)
speed_all = load([PathName, FileName]);
num_ds = speed_all(end,1);



%calculate histogram global switch ("gs")
speed_gs = zeros(size(speed_all,1),4);
speed_high = zeros(size(speed_all,1),4);
speed_low = zeros(size(speed_all,1),4);
num_tr_gs = 0;
num_tr_high = 0;
num_tr_low = 0;

for i=1:num_ds
    ind_s = speed_all(:,1)==i;
    speed_s = speed_all(ind_s,2:4);
    last_tr = speed_s(end,1);
    h = waitbar(0,['Dataset: ', int2str(i),'. Calculate velocities ...']);
    for j=1:last_tr
        tr_id = find(speed_s(:,1)==j);
        if ~isempty(tr_id)
            time_tr = speed_s(tr_id, 2);
            ind_t = find(time_tr>min_t & time_tr<max_t);
            if ~isempty(ind_t)
                num_tr_gs = num_tr_gs + 1;
                dp = length(ind_t);
                null_ind = find(speed_gs==0,1);
                ind = null_ind:null_ind+dp-1;
                speed_gs(ind,1) = i;
                speed_gs(ind,2:4) = speed_s(ind_t + tr_id(1) - 1, :);
            end
            ind_t = find(time_tr < min_t);
            if ~isempty(ind_t)
                num_tr_high = num_tr_high + 1;
                dp = length(ind_t);
                null_ind = find(speed_high(:,1)==0,1);
                ind = null_ind:null_ind+dp-1;
                speed_high(ind,1) = i;
                speed_high(ind,2:4) = speed_s(ind_t + tr_id(1) - 1,:);
            end
            ind_t = find(time_tr > max_t);
            if ~isempty(ind_t)
                num_tr_low = num_tr_low + 1;
                dp = length(ind_t);
                null_ind = find(speed_low(:,1)==0,1);
                ind = null_ind:null_ind+dp-1;
                speed_low(ind,1) = i;
                speed_low(ind,2:4) = speed_s(ind_t + tr_id(1) - 1,:);
            end
            waitbar(j/last_tr);
        end
    end
    close(h)
end

speed_gs(find(speed_gs(:,1)==0,1):end,:)=[];
speed_high(find(speed_high(:,1)==0,1):end,:)=[];
speed_low(find(speed_low(:,1)==0,1):end,:)=[];


if ~isempty(speed_gs)
    assignin('base', 'speed_gs', speed_gs);
    assignin('base', 'num_tr_gs', num_tr_gs);
    path_gs = [PathName, 'Speed_gs.txt'];
    save(path_gs , 'speed_gs', '-ascii', '-double' ,'-tabs'); 
    path_gs = [PathName, 'num_tr_gs.txt'];
    save(path_gs , 'num_tr_gs', '-ascii', '-double' ,'-tabs'); 
    figure(1), 
    set(1, 'Name', 'Global Switching Window');
    [n,xout] = hist(speed_gs(:,4), 0:1/20:5);
    subplot(2,1,1),
    bar(xout,n, .8, 'r'), 
    %fitting a gaussian to speed histogram
    f = fittype('gauss2');
    fit_start = 3;
    res_fit  = fit(xout(fit_start:end)',n(fit_start:end)',f);
    confi_speed = confint(res_fit);
    err_speed_1 = confi_speed(2,2) - confi_speed(1,2);
    err_speed_2 = confi_speed(2,5) - confi_speed(1,5);
    title('Speed Histogram' , 'FontSize', 14),
    hold on
    plot(res_fit, '-k'),
    hold off
    xlabel('Speed [µm/s]', 'FontSize', 12), ...
    ylabel('Frequency', 'FontSize', 12),
    xlim([0 3]);
    grid on;
    text(2.5,max(n), ['\Rightarrow v_1=(',...
        num2str(res_fit.b1,'%1.3f'),'\pm', num2str(err_speed_1, ...
        '%1.3f)µm/s'),'\newline \Rightarrow v_2=(',...
        num2str(res_fit.b2,'%1.3f'),'\pm', num2str(err_speed_2, ...
        '%1.3f)µm/s')], ...
         'FontSize',12),
    subplot(2,1,2),
    time = speed_gs(:,3);
    plot(time , speed_gs(:,4), '+r');
    title('Speed versus Time', 'FontSize', 14),
    xlabel('Time [min]', 'FontSize', 12), ...
    ylabel('Speed [µm/s]', 'FontSize', 12),
    xlim([min(time) max(time)]);
    ylim([0 3]);
    grid on;
end


if ~isempty(speed_high)
    assignin('base', 'speed_high', speed_high);
    assignin('base', 'num_tr_high', num_tr_high);
    path_high = [PathName, 'Speed_high.txt'];
    save(path_high , 'speed_high', '-ascii', '-double' ,'-tabs'); 
    path_high = [PathName, 'num_tr_high.txt'];
    save(path_high , 'num_tr_high', '-ascii', '-double' ,'-tabs'); 
    figure(2),    
    set(2, 'Name', 'High Velocity Mode');
    [n,xout] = hist(speed_high(:,4), 0:1/20:5);
    subplot(2,1,1),
    bar(xout,n, .8, 'r'), 
    %fitting a gaussian to speed histogram
    f = fittype('gauss1');
    fit_start = 8;
    res_fit  = fit(xout(fit_start:end)',n(fit_start:end)',f);
    confi_speed = confint(res_fit);
    err_speed = confi_speed(2,2) - confi_speed(1,2);
    title('Speed Histogram' , 'FontSize', 14),
    hold on
    plot(res_fit, '-k'),
    hold off
    xlabel('Speed [µm/s]', 'FontSize', 12), ...
    ylabel('Frequency', 'FontSize', 12),
    xlim([0 3]);
    grid on;
    text(1.8,max(n), ['\Rightarrow v =(',num2str(res_fit.b1,'%1.3f'),...
        '\pm', num2str(err_speed, '%1.3f)µm/s'), '\newline Mean velocity =', ...
        num2str(mean(speed_high(:,4)), '%1.3fµm/s')], 'FontSize',12),
    subplot(2,1,2),
    time = speed_high(:,3);
    plot(time , speed_high(:,4), '+r');
    title('Speed versus Time', 'FontSize', 14),
    xlabel('Time [min]', 'FontSize', 12), ...
    ylabel('Speed [µm/s]', 'FontSize', 12),
    xlim([min(time) max(time)]);
    ylim([0 3]);
    grid on;
end

if ~isempty(speed_low)
    assignin('base', 'speed_low', speed_low);
    assignin('base', 'num_tr_low', num_tr_low);
    path_low = [PathName, 'Speed_low.txt'];
    save(path_low , 'speed_low', '-ascii', '-double' ,'-tabs'); 
    path_low = [PathName, 'num_tr_low.txt'];
    save(path_low , 'num_tr_low', '-ascii', '-double' ,'-tabs'); 
    figure(3),    
    set(3, 'Name', 'Low Velocity Mode');
    [n,xout] = hist(speed_low(:,4), 0:1/20:5);
    subplot(2,1,1),
    bar(xout,n, .8, 'r'), 
    %fitting a gaussian to speed histogram
    f = fittype('gauss1');
    fit_start = 5;
    res_fit  = fit(xout(fit_start:end)',n(fit_start:end)',f);
    confi_speed = confint(res_fit);
    err_speed = confi_speed(2,2) - confi_speed(1,2);
    title('Speed Histogram' , 'FontSize', 14),
    hold on
    plot(res_fit, '-k'),
    hold off
    xlabel('Speed [µm/s]', 'FontSize', 12), ...
    ylabel('Frequency', 'FontSize', 12),
    xlim([0 3]);
    grid on;
    text(1.8,max(n), ['\Rightarrow v =(',num2str(res_fit.b1,'%1.3f'),...
        '\pm', num2str(err_speed, '%1.3f)µm/s'), '\newline Mean velocity =', ...
        num2str(mean(speed_low(:,4)), '%1.3fµm/s')], 'FontSize',12),
    subplot(2,1,2),
    time = speed_low(:,3);
    plot(time , speed_low(:,4), '+r');
    title('Speed versus Time', 'FontSize', 14),
    xlabel('Time [min]', 'FontSize', 12), ...
    ylabel('Speed [µm/s]', 'FontSize', 12),
    xlim([min(time) max(time)]);
    ylim([0 3]);
    grid on;
end
