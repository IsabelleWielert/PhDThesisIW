function reanalyze_data(gaussian)

%set graphical objects for histogram
h_fig = figure(1);
set(h_fig, 'Name', 'Speed histogram');
h_ax = axes;
h_hist = bar(h_ax,0,0,.8,'r');
hist_XLabel = get(h_ax, 'XLabel');
set(hist_XLabel, 'String', 'Speed [µm/s]');
hist_YLabel = get(h_ax, 'YLabel');
set(hist_YLabel, 'String', 'Frequency');
hist_Title = get(h_ax, 'Title');
set(hist_Title, 'String', 'Speed histogram', 'FontWeight', 'bold');
set(h_ax, 'XLim', [0 3]);
grid(h_ax, 'on');

%open dialog for choosing first data_all.txt
[FileName,PathName] = uigetfile('*data_all.txt');
%if cancel is pushed, then return.
if (isscalar(FileName) == 1) && (isscalar(PathName) == 1);
    return;
end  

cd(PathName)
data_paths = dir('*tr_all.txt');
num_tr_paths = dir('*num_tracks.txt');
num_files = size(data_paths,1);

data_all = zeros(1000000,6);
dataset_id = 0;
num_tr = 0;
for i=1:num_files
    data = load(data_paths(i).name);
    num_tr = num_tr + load(num_tr_paths(i).name);
    null_id = find(data_all(:,1)==0,1);
    ind = null_id:null_id + size(data,1) - 1;
    data_all(ind,1) = data(:,1) + dataset_id;
    data_all(ind,2:6) = data(:,2:6);
    dataset_id = max(data_all(:,1));
end
data_all(find(data_all(:,1)==0,1):end,:)=[];


%calculate speed histogram
x_hist = 0:1/20:5;
speed = data_all(:,6);
assignin('base', 'speed', speed)
n = hist(speed, x_hist);
n = n./sum(n);
%plot speed histogram
set(h_hist, 'XData', x_hist, 'YData', n); 
speed_av = mean(speed);
speed_av_err = std(speed)/sqrt(num_tr);

%fitting a gaussian to speed histogram

if (gaussian == 1)
    f = fittype('gauss1');
    start_ind = 25;
elseif (gaussian == 2)
    f = fittype('gauss2');
    start_ind = 2;
end
res_fit  = fit(x_hist(start_ind:end)',n(start_ind:end)',f);
confi_speed = confint(res_fit);

if (gaussian == 1)
     err_speed = confi_speed(2,2) - confi_speed(1,2);
elseif (gaussian == 2)
     err_speed_1 = confi_speed(2,2) - confi_speed(1,2);
     err_speed_2 = confi_speed(2,5) - confi_speed(1,5);
end

hold(h_ax, 'on')
plot(h_ax, x_hist, res_fit(x_hist), '-k'),
hold(h_ax, 'off')

y_lim = get(h_ax, 'YLim');
y_max = y_lim(end);
if (gaussian == 1)
    text(.1,y_max, ...
    ['Single Gaussian Fit: \newline\Rightarrow v =(',...
        num2str(res_fit.b1,'%1.3f'),'\pm', num2str(err_speed, ...
        '%1.3f)µm/s'), '\newline Arithmetic Mean: \newline\Rightarrow v =(', ...
        num2str(speed_av, '%1.3f'),'\pm', num2str(speed_av_err, ...
        '%1.3f)µm/s'), '\newline# of tracks:', int2str(num_tr)],...
        'VerticalAlignment', 'top'),
elseif (gaussian == 2)
text(.1,y_max, ...
['Double Gaussian Fit\newline\Rightarrow v_1=(',...
    num2str(res_fit.b1,'%1.3f'),'\pm', num2str(err_speed_1, ...
    '%1.3f)µm/s'),'\newline\Rightarrow v_2=(',...
    num2str(res_fit.b2,'%1.3f'),'\pm', num2str(err_speed_2, ...
    '%1.3f)µm/s'),'\newline# of tracks:', int2str(num_tr)],...
    'VerticalAlignment', 'top'),
end

assignin('base', 'data_all', data_all)
[FileName,PathName] = uiputfile('*.txt');
%if cancel is pushed, then return.
if (isscalar(FileName) == 1) && (isscalar(PathName) == 1);
    return;
end  
path_data = [PathName, FileName]; 
path_num_tr = [PathName, FileName(1:end-4), '_num_tracks.txt']; 
save(path_data , 'data_all', '-ascii', '-double' ,'-tabs'); 
save(path_num_tr , 'num_tr', '-ascii'); 


function speed = calc_speed(pos, time, dist)

tr_len = length(time);
dp = tr_len - dist;
speed = zeros(dp,1);
for j=1:dp;
    ind = j + dist;
    if (ind <= tr_len)
        s_vec = [pos(ind,1)-pos(j,1);pos(ind,2)-pos(j,2)];
        delta_t = time(ind)-time(j);   
        speed(j) = norm(s_vec)/delta_t;
    end
end