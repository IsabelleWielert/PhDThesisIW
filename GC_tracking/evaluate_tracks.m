 function varargout = evaluate_tracks(varargin)
% EVALUATE_TRACKS M-file for evaluate_tracks.fig
%      EVALUATE_TRACKS, by itself, creates a new EVALUATE_TRACKS or raises the existing
%      singleton*.
%
%      H = EVALUATE_TRACKS returns the handle to a new EVALUATE_TRACKS or the handle to
%      the existing singleton*.
%
%      EVALUATE_TRACKS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EVALUATE_TRACKS.M with the given input arguments.
%
%      EVALUATE_TRACKS('Property','Value',...) creates a new EVALUATE_TRACKS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before evaluate_tracks_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to evaluate_tracks_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help evaluate_tracks

% Last Modified by GUIDE v2.5 07-Nov-2011 11:04:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @evaluate_tracks_OpeningFcn, ...
                   'gui_OutputFcn',  @evaluate_tracks_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before evaluate_tracks is made visible.
function evaluate_tracks_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to evaluate_tracks (see VARARGIN)

% Choose default command line output for evaluate_tracks
handles.output = hObject;

camera = varargin{1};
binning = varargin{3};
% pixel factor for absolute lengthscales
switch camera
    case 1 %pixel scaling in units of micron for 
           %100x oil immersion objective (TIRF setup)
           handles.magni=binning*0.083;
        set(handles.cam_name, 'String', 'TIRF');
    case 2 %pixel scaling in units of micron for 
           %100x oil immersion objective (Confocal setup)
        handles.magni = binning*0.065;
        set(handles.cam_name, 'String', 'Confocal');
    case 3  %pixel scaling in units of micron for 
            %100x oil immersion objective (Zeiss OT setup)
        handles.magni = binning*0.099;
        set(handles.cam_name, 'String', 'Zeiss');
    case 4  %pixel scaling in units of micron for 
            %100x oil immersion objective (Attic lab OT setup)
        handles.magni = binning*0.059;
        set(handles.cam_name, 'String', 'Zeiss');
end    
set(handles.scale, 'String', num2str(handles.magni, '%1.3f µm/Pixel'));
handles.tif_format = varargin{2};
handles.fps = varargin{4};
set(handles.fps_out, 'String', int2str(handles.fps));
handles.delta_t = varargin{5};
set(handles.delta_t_in, 'String', num2str(handles.delta_t, '%1.1f'));
handles.dist = round(handles.delta_t*handles.fps);

[FileName_tr,PathName_tr] = uigetfile('*.txt', ...
                                    'Choose txt-file containing tracks');
if (isscalar(FileName_tr) == 1) && (isscalar(PathName_tr) == 1);
    return;
end 

cd(PathName_tr)
handles.tr = load(FileName_tr);

cd('..')
[FileName_tif,PathName_tif] = uigetfile('*.tif', ...
                             'Choose Tif-file of corresponding movie');
if (isscalar(FileName_tif) == 1) && (isscalar(PathName_tif) == 1);
    guidata(hObject, handles);
    return;
end 
handles.img_path = {PathName_tif, FileName_tif};
cd(PathName_tif);
drawnow

if (handles.tif_format <= 2)
    handles.single_tifs = dir('*.tif');
    handles.frames = size(handles.single_tifs,1);
    handles.tif_info = imfinfo([handles.img_path{1}, handles.img_path{2}]);
    handles.img_width = handles.tif_info(1).Width; %(1) hinzugefügt
    handles.img_height = handles.tif_info(1).Height; %(1) hinzugefügt
else
    handles.tif_info = imfinfo([handles.img_path{1}, handles.img_path{2}]);
    handles.frames = size(handles.tif_info,1);
    handles.img_width=handles.tif_info(1).Width;
    handles.img_height=handles.tif_info(1).Height;
end
set(handles.frames_out, 'String', int2str(handles.frames));
set(handles.time_out, 'String', num2str(handles.frames/handles.fps, ...
                                                                 '%3.1f'));
 
cd(PathName_tr);
handles.file_path = cd;

handles.num_tracks = max(handles.tr(:,4));
set(handles.num_tracks_out, 'String', int2str(handles.num_tracks));
%calc cell density
img_area = handles.img_width*handles.img_height*(handles.magni^2);
cell_dens = calc_cell_density(handles.tr, 1:max(handles.tr(:,3)), img_area);
set(handles.total_cell_dens_out, 'String', num2str(cell_dens, '%1.4f'));
drawnow

%for msd calculations 
handles.my_msd = fittype('2.*tau_c*(v_c^2).*(tau-tau_c.*(1-exp(-tau./tau_c)))+A',...
    'independent', 'tau');
A0 = 0.05;
tau_c0 = 1.5;
v_c0 = 1.5;
handles.opts_msd = fitoptions(handles.my_msd);
set(handles.opts_msd,'TolFun',1E-6, 'TolX', 1E-6, 'StartPoint', ...
    [A0, tau_c0, v_c0]); 


%for acf calculations  
handles.my_acf = fittype('v_c^2*exp(-tau./tau_c)', 'independent', 'tau');
%intial values for non linear fit
tau_c0 = 1.5;
v_c0 = 1.5; 
handles.opts_acf = fitoptions(handles.my_acf);
set(handles.opts_acf,'TolFun',1E-6, 'TolX', 1E-6, 'StartPoint', ...
    [tau_c0, v_c0]); 

%set graphical objects for histogram
set(0, 'CurrentFigure', handles.evaluation);
movegui(handles.evaluation, 'northwest')
set(handles.evaluation, 'CurrentAxes', handles.hist_axes)
handles.hist = bar(handles.hist_axes,0,0,.8,'r');
hist_XLabel = get(handles.hist_axes, 'XLabel');
set(hist_XLabel, 'String', 'Speed [µm/s]');
hist_YLabel = get(handles.hist_axes, 'YLabel');
set(hist_YLabel, 'String', 'Frequency');
hist_Title = get(handles.hist_axes, 'Title');
set(hist_Title, 'String', 'Speed histogram', 'FontWeight', 'bold');
set(handles.hist_axes, 'XLim', [0 3]);
grid(handles.hist_axes, 'on');

%set graphical objects for speed vs time
set(handles.evaluation, 'CurrentAxes', handles.speed_vs_time_axes)
handles.speed_time = plot(handles.speed_vs_time_axes, 0, 0, '-r');
hold(handles.speed_vs_time_axes, 'on');
handles.time_ind = plot(handles.speed_vs_time_axes, [0,0], [0,3], '-k');
hold(handles.speed_vs_time_axes, 'off');
hist_XLabel = get(handles.speed_vs_time_axes, 'XLabel');
set(hist_XLabel, 'String', 'Time [s]');
hist_YLabel = get(handles.speed_vs_time_axes, 'YLabel');
set(hist_YLabel, 'String', 'Speed [µm/s]');
hist_Title = get(handles.speed_vs_time_axes, 'Title');
set(hist_Title, 'String', 'Speed versus Time','FontWeight', 'bold');
set(handles.speed_vs_time_axes, 'XLim', [0 handles.frames/handles.fps],...
                                                             'YLim',[0,3]);
grid(handles.speed_vs_time_axes, 'on')

%set graphical objects for movies
set(handles.evaluation, 'CurrentAxes', handles.movie_axes)
colormap(handles.movie_axes, 'gray');
set(handles.movie_axes, 'CLim', [0 200]);
img = zeros(handles.img_height, handles.img_width);
handles.img_mov = image(img, 'CDataMapping', 'scaled');
handles.rect = rectangle('Position',[0,0,handles.img_width, handles.img_height],...
        'LineWidth', 2, 'LineStyle', '--');
axis(handles.movie_axes, 'image');
hold(handles.movie_axes, 'on')
handles.track_mov = plot(handles.movie_axes, 1, 1, '-k');
hold(handles.movie_axes, 'on')
handles.cell_mov = plot(handles.movie_axes, handles.img_width/2, handles.img_height/2,...
      'o', 'MarkerEdgeColor','k', 'MarkerFaceColor','r',  'MarkerSize', 3);
hold(handles.movie_axes, 'off');  
                    
set(handles.evaluation, 'CurrentAxes', handles.track_movie_axes)
colormap(handles.track_movie_axes, 'gray');
set(handles.track_movie_axes, 'CLim', [0 200]);
set(handles.zoom_rg, 'String', '120');
handles.img_zoom =120; 
img = zeros(handles.img_zoom, handles.img_zoom);
handles.img_tr_mov = image(img, 'CDataMapping', 'scaled');
axis(handles.track_movie_axes, 'image');
hold(handles.track_movie_axes, 'on')
handles.track_tr_mov = plot(handles.track_movie_axes, 1, 1, '-k');
hold(handles.movie_axes, 'on')
handles.cell_tr_mov = plot(handles.track_movie_axes, handles.img_zoom/2,...
   handles.img_zoom/2, 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',...
    'r',  'MarkerSize', 3);
hold(handles.track_movie_axes, 'off');

%define figure for msd output
handles.msd_plot_fig = figure;
movegui(handles.msd_plot_fig, 'northeast')
set(handles.msd_plot_fig,'NumberTitle', 'off', 'Name', 'Mean Square Displacement');
set(handles.msd_plot_fig,'Visible', 'off');
handles.msd_plot_axes = axes;

handles.msd_plot = errorbar(handles.msd_plot_axes, 0, 0, 0, 0, '+k');
hold(handles.msd_plot_axes, 'on')
handles.msd_fit_plot = plot(handles.msd_plot_axes,0,0,'-r');
hold(handles.msd_plot_axes, 'off');
grid(handles.msd_plot_axes, 'on');

msd_XLabel = get(handles.msd_plot_axes, 'XLabel');
set(msd_XLabel, 'String', 'time \tau [s]');
msd_YLabel = get(handles.msd_plot_axes, 'YLabel');
set(msd_YLabel, 'String', 'MSD [µm^2]');
msd_Title = get(handles.msd_plot_axes, 'Title');
set(msd_Title, 'String', 'Mean Square Displacement (MSD)', 'FontWeight', 'bold');

%define figure for acf output
handles.acf_plot_fig = figure;
movegui(handles.acf_plot_fig, 'southeast')
set(handles.acf_plot_fig,'NumberTitle', 'off', 'Name', 'Angle Correlation Function');
set(handles.acf_plot_fig,'Visible', 'off');
handles.acf_plot_axes = axes;

handles.acf_plot = errorbar(handles.acf_plot_axes, 0, 0, 0, 0, '+k');
hold(handles.acf_plot_axes, 'on')
handles.acf_fit_plot = plot(handles.acf_plot_axes, 0, 0, '-r');
hold(handles.acf_plot_axes, 'off');
grid(handles.acf_plot_axes, 'on');

acf_XLabel = get(handles.acf_plot_axes, 'XLabel');
set(acf_XLabel, 'String', 'time \tau [s]');
acf_YLabel = get(handles.acf_plot_axes, 'YLabel');
set(acf_YLabel, 'String', 'ACF');
acf_Title = get(handles.acf_plot_axes, 'Title');
set(acf_Title, 'String', 'Angle Correlation Function (ACF)', 'FontWeight', 'bold');

handles.tr_id = 1;
handles.total_dp = 0;
handles.num_dead = 0;
set(handles.curr_tr, 'String', int2str(handles.tr_id));
drawnow
guidata(hObject, handles);
evaluate_track(handles);

 

% UIWAIT makes evaluate_tracks wait for user response (see UIRESUME)
% uiwait(handles.evaluation);


% --- Outputs from this function are returned to the command line.
function varargout = evaluate_tracks_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




function evaluate_track(handles)

assignin('base', 'h', handles.evaluation)
%get position data for single track
ind_tr = find(handles.tr(:,4)==handles.tr_id);
pos = handles.tr(ind_tr,1:2);

%get frame identifier
fr_id = handles.tr(ind_tr,3);

%calc cell density for track
%img_area in units of mm^2
img_area = handles.img_width*handles.img_height*(handles.magni)^2;
cell_dens_tr = calc_cell_density(handles.tr, fr_id, img_area);
set(handles.cell_dens_out, 'String', num2str(cell_dens_tr, '%1.4f'));

%fill pos data of dead frames via linear interpolation
pos = fill_dead_frames(fr_id, pos);

%reset fr_id to single frame steps
fr_id = fr_id(1):fr_id(end);
dp = length(fr_id) - handles.dist;
time = 1/handles.fps*(fr_id-1); %units of seconds

%calc speed via pos data
speed = calc_speed(pos, time, handles.dist);

%calculate correlation time via msd for time range defined by min_tau and
%max_tau (mininal track lenght: 1min)
if dp/handles.fps >= 20 %vorher stand hier 60 isa 
    tau_max = 15; %in units of seconds
    if get(handles.plot_msd, 'Value')
        if tau_max > time(end-handles.dist-1)-time(1);
            tau_max = time(end-handles.dist-1)-time(1);
        end

        tau_fr = 1:round(handles.fps*tau_max); %tau in units of frames
        tau = tau_fr/handles.fps; %tau in units of seconds
        [msd, msd_err] = calc_msd(pos, tau_fr(end));
        set(handles.opts_msd, 'Weights', 1./(msd_err).^2);

        if (sum(isnan(msd))==0)
            [msd_fun, msd_gof] = fit(tau', msd, handles.my_msd, handles.opts_msd);
            tau_c = msd_fun.tau_c;
            msd_rsquare = msd_gof.rsquare;
            set(handles.msd_tau_c_out, 'String', num2str(tau_c, '%2.2f'));
            set(handles.msd_gof_out, 'String', num2str(msd_rsquare, '%1.4f'));
            set(handles.msd_plot, 'XData', tau, 'YData', msd, 'LData', msd_err, ...
                'UData', msd_err),
            set(handles.msd_fit_plot , 'XData', tau, 'YData', msd_fun(tau));
        else
            set(handles.msd_tau_c_out, 'String', 'NaN');
            set(handles.msd_gof_out, 'String', 'NaN');
        end
    else
        set(handles.msd_tau_c_out, 'String', 'NaN');
        set(handles.msd_gof_out, 'String', 'NaN');
    end

    %calc correlation time of directed movement via velocity auto correlation 
    % function
    if get(handles.plot_acf, 'Value')
        [acf, acf_err] = calc_acf(pos, time, tau_fr(end), handles.dist);
        set(handles.opts_acf, 'Weights', 1./(acf_err).^2);
        %do exponetial fit with my_acf;
        if (sum(isnan(acf))==0)
            [acf_fun, acf_gof] = fit(tau', acf, handles.my_acf, handles.opts_acf); 
            tau_c = acf_fun.tau_c;
            acf_rsquare = acf_gof.rsquare;
            set(handles.acf_tau_c_out, 'String', num2str(tau_c, '%2.2f'));
            set(handles.acf_gof_out, 'String', ...
                                                   num2str(acf_rsquare, '%1.4f'));
            set(handles.acf_plot, 'XData', tau, 'YData', acf, ...
                'LData', acf_err, 'UData', acf_err),
            set(handles.acf_fit_plot , 'XData', tau, 'YData', acf_fun(tau));
        else
            set(handles.acf_tau_c_out, 'String', 'NaN');
            set(handles.acf_gof_out, 'String', 'NaN');
        end
    else
        set(handles.acf_tau_c_out, 'String', 'NaN');
        set(handles.acf_gof_out, 'String', 'NaN');
    end
else
    set(handles.msd_tau_c_out, 'String', 'NaN');
    set(handles.msd_gof_out, 'String', 'NaN');
    set(handles.acf_tau_c_out, 'String', 'NaN');
    set(handles.acf_gof_out, 'String', 'NaN');
end
                                            
%set outputs in gui
set(handles.first_fr_out, 'String', int2str(fr_id(1)));
handles.mov_1st_frame = 1;
set(handles.mov_1st_frame_out, 'String', '1');
set(handles.frame_slider, 'Min', 1,'Max', dp,...
    'Value', 1, 'SliderStep', [1/(dp-1), 10/(dp-1)]);
set(handles.last_fr_out, 'String', int2str(fr_id(end)));
set(handles.time_tr_out, 'String', num2str(time(end)-time(1), '%3.1f'));
set(handles.dp_out, 'String', int2str(dp));
set(handles.total_dp_out, 'String', int2str(handles.total_dp));
drawnow

    
%calculate speed histogram
x_hist = 0:1/20:5;
n = hist(speed, x_hist);

%remove last "dist" elements of time vector and pos list to adjust 
%length to speed vector
time(dp+1:end)=[];
pos(dp+1:end,:)=[];

%rescale pos list in units of pixel
pos_pxl = pos/handles.magni;

%plot speed histogram
set(handles.hist, 'XData', x_hist, 'YData', n./sum(n)); 

%plot speed versus time
set(handles.speed_time, 'XData', time, 'YData', speed);
set(handles.time_ind, 'XData', [time(1), time(1)]);

%load first frame of movie
switch handles.tif_format
    case 1
        img = imread([handles.img_path{1},...
                        handles.single_tifs(fr_id(1)).name]); 
    case 2
        img = imread([handles.img_path{1},...
                        handles.single_tifs(fr_id(1)).name]); 
        %conversion to 8-bit
        img = uint8(round(img./(2^8)));
    case 3
        img = imread([handles.img_path{1}, handles.img_path{2}], ...
                          fr_id(1));%, 'Info', handles.tif_info);
    case 4  
        img = imread([handles.img_path{1}, handles.img_path{2}], ...
                          fr_id(1));%, 'Info', handles.tif_info);
        %conversion to 8-bit
        img = uint8(round(img./(2^8)));
end  


%show first img in axes "movie_axes"
%range of track in units of pixel
x_min = floor(min(pos_pxl(:,1)))-10;
y_min = floor(min(pos_pxl(:,2)))-10;
x_max = ceil(max(pos_pxl(:,1)))+10;
y_max = ceil(max(pos_pxl(:,2)))+10;
set(handles.img_mov, 'CData', img)
set(handles.rect, 'Position',[x_min,y_min,x_max-x_min,y_max-y_min]);
set(handles.cell_mov, 'XData', pos_pxl(1,1), 'YData', pos_pxl(1,2));
set(handles.track_mov, 'XData', pos_pxl(:,1), 'YData', pos_pxl(:,2));                                         

%show first img i n axes "movie_axes"
x_min_tr = floor(min(pos_pxl(1,1))) - round(handles.img_zoom/2);
y_min_tr = floor(min(pos_pxl(1,2))) - round(handles.img_zoom/2);
x_max_tr = ceil(max(pos_pxl(1,1))) + round(handles.img_zoom/2);
y_max_tr = ceil(max(pos_pxl(1,2))) + round(handles.img_zoom/2);
if x_min_tr <= 1;
x_min_tr = 1;
end
if x_max_tr >= handles.img_width;
    x_max_tr = handles.img_width;
end
if y_min_tr <= 1;
    y_min_tr = 1;
end
if y_max_tr >= handles.img_height
    y_max_tr = handles.img_height;
end
set(handles.img_tr_mov, 'CData', img(y_min_tr:y_max_tr, ...
                                                   x_min_tr:x_max_tr));
set(handles.cell_tr_mov, 'XData', pos_pxl(1,1) - x_min_tr + 1, 'YData',...
                                           pos_pxl(1,2)- y_min_tr + 1);
set(handles.track_tr_mov, 'XData', pos_pxl(:,1)- x_min_tr + 1, 'YData',...
                                           pos_pxl(:,2)- y_min_tr + 1);
set(handles.track_movie_axes, 'XLim', [1, x_max_tr-x_min_tr], ...
                                       'YLim', [1, y_max_tr-y_min_tr]); 
set(handles.frame_id, 'String', num2str(fr_id(1), '%4.0f'));                                     
set(handles.time_id, 'String', num2str(time(1), '%3.1fs'));                                   
%update handles
handles.dp = dp;
handles.data = zeros(dp, 5);
handles.pos = pos;
handles.time = time;
handles.speed = speed;
handles.fr_id = fr_id;
guidata(handles.evaluation, handles); 

% first_speed_check calculates a moving average of the speed vector 
% to determine motility pauses or too high velocities
%red Background in curr_tr stands for bad track
set(handles.curr_tr, 'BackgroundColor', 'white')
win = handles.fps*5; %in units of frames
min_speed = 0.15; %in units of µm/s
max_speed = 5; %in units of µm/s
bad_tr = first_speed_check(speed, win, min_speed, max_speed);
if bad_tr == 1
    set(handles.curr_tr, 'BackgroundColor', 'red') 
elseif  get(handles.auto, 'Value') == 1
    save_Callback(handles.save,0, handles)
end



function pos = fill_dead_frames(fr_id, pos)

dp = length(fr_id);
pos_old = pos;
num_filled = 0;
for i=2:dp
    fr_diff = fr_id(i)-fr_id(i-1);
    if fr_diff>1
        pos_diff = pos_old(i,:) - pos_old(i-1,:);
        pos_filled = zeros(fr_diff-1,2);
        for j=1:fr_diff-1
            pos_filled(j,:) = pos_old(i-1,:)+j/fr_diff.*pos_diff;         
        end
        ind_fill = num_filled+i-1;
        pos_tmp1 = pos(1:ind_fill,:);
        pos_tmp2 = pos(ind_fill+1:end,:);
        pos = vertcat(pos_tmp1, pos_filled, pos_tmp2);
        num_filled = num_filled + fr_diff-1;
    end
end

function cell_dens = calc_cell_density(tr, fr_id, img_area)

dp = length(fr_id);
num_cells = zeros(dp,1);
for i=1:dp
    num_cells(i) = length(find(tr(:,3)==fr_id(i)));
end
cell_dens = mean(num_cells)/img_area; 


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

function [acf, acf_err] = calc_acf(pos, time, max_tau, dist)

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

acf = zeros(max_tau,1);
acf_err = zeros(max_tau,1);
for j=1:max_tau;
    acf_j = zeros(dp-j,1);
    for l=1:dp-j
        ind = l+j;
        acf_j(l) = (v(l,1)*v(ind,1)+v(l,2)*v(ind,2));
    end
    if ~isempty(acf_j)
        acf(j) = mean(acf_j);
        acf_err(j) = std(acf_j)/sqrt(dp-j);
    end
end

function [msd, msd_err] = calc_msd(pos,max_tau)

tr_len = size(pos,1);
msd = zeros(max_tau,1);
msd_err = zeros(max_tau,1);

for j=1:max_tau;
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
    end
end 

function bad_tr = first_speed_check(speed, win, min_speed, max_speed)

bad_tr = 0;
dp = length(speed);

if dp >= win
      for k=1:dp-win
         speed_av = mean(speed(k:k+win));
         if (speed_av < min_speed || speed_av > max_speed) 
             bad_tr = 1;
             break;
         end
      end
end






% --- Executes on button press in next_tr.
function next_tr_Callback(hObject, eventdata, handles)
% hObject    handle to next_tr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.tr_id < handles.num_tracks)
    handles.tr_id = handles.tr_id + 1;
end
set(handles.curr_tr, 'String', int2str(handles.tr_id));
guidata(hObject, handles);
evaluate_track(handles);

% --- Executes on button press in play_tr.
function play_tr_Callback(hObject, eventdata, handles)
% hObject    handle to play_tr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of play_tr

while(get(hObject,'Value')==1)
    i = handles.mov_1st_frame;
    pos_pxl = handles.pos/handles.magni;
    skip_frames = str2double(get(handles.skip_frames, 'String'));    
    while(get(hObject,'Value')==1 && i<= handles.dp)
        update_movie(pos_pxl, handles, i)
        i = i + skip_frames;
    end
    if i>handles.dp
        set(hObject, 'Value', 0)
        handles.mov_1st_frame = 1;
        update_movie(pos_pxl, handles, 1)
    else
        handles.mov_1st_frame = i-skip_frames;
    end
    set(handles.frame_slider, 'Value', handles.mov_1st_frame);
    set(handles.mov_1st_frame_out, 'String', int2str(handles.mov_1st_frame));
    guidata(handles.evaluation, handles);
end


function update_movie(pos_pxl, handles, i)

switch handles.tif_format
   case 1
        img = imread([handles.img_path{1},...
                        handles.single_tifs(handles.fr_id(i)).name]); 
    case 2
        img = imread([handles.img_path{1},...
                        handles.single_tifs(handles.fr_id(i)).name]); 
        %conversion to 8-bit
        img = uint8(round(img./(2^8)));
    case 3
        img = imread([handles.img_path{1}, handles.img_path{2}], ...
                          handles.fr_id(i));%, 'Info', handles.tif_info);
    case 4  
        img = imread([handles.img_path{1}, handles.img_path{2}], ...
                          handles.fr_id(i));%, 'Info', handles.tif_info);
        %conversion to 8-bit
        img = uint8(round(img./(2^8)));
end  
%show img in axes "movie_axes"
%     set(handles.img_mov, 'CData', img)
set(handles.cell_mov, 'XData', pos_pxl(i,1), 'YData', pos_pxl(i,2));

%show img in axes "movie_axes"
x_min_tr = floor(min(pos_pxl(i,1)) - handles.img_zoom/2);
y_min_tr = floor(min(pos_pxl(i,2)) - handles.img_zoom/2);
x_max_tr = ceil(max(pos_pxl(i,1)) + handles.img_zoom/2);
y_max_tr = ceil(max(pos_pxl(i,2)) + handles.img_zoom/2);
if x_min_tr <= 1;
x_min_tr = 1;
end
if x_max_tr >= handles.img_width;
    x_max_tr = handles.img_width;
end
if y_min_tr <= 1;
    y_min_tr = 1;
end
if y_max_tr >= handles.img_height
    y_max_tr = handles.img_height;
end
set(handles.img_tr_mov, 'CData', img(y_min_tr:y_max_tr, ...
                                                   x_min_tr:x_max_tr));
set(handles.cell_tr_mov, 'XData', pos_pxl(i,1) - x_min_tr + 1, 'YData',...
                                           pos_pxl(i,2)- y_min_tr + 1);
set(handles.track_tr_mov, 'XData', pos_pxl(:,1)- x_min_tr + 1, 'YData',...
                                           pos_pxl(:,2)- y_min_tr + 1);
set(handles.track_movie_axes, 'XLim', [1, x_max_tr-x_min_tr], ...
                                       'YLim', [1, y_max_tr-y_min_tr]);
%update time indicator in speed_vs_time_axes and in time_id edit
set(handles.time_ind, 'XData', [handles.time(i),handles.time(i)]);
set(handles.frame_id, 'String', num2str(handles.fr_id(i), '%4.0f'));     
set(handles.time_id, 'String', num2str(handles.time(i), '%3.1fs'));
drawnow
% frame = getframe(handles.track_movie_axes);
% [img,~] = frame2im(frame);
% %img = rgb2gray(img);
% %generation of a 8-bit multi-image tiff-file with tracks
% imwrite(img, 'GC_plus_track.tif', 'tif', ...
%         'WriteMode', 'append', 'Compression', 'none')
%pause(0.001);

% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%store graphical outputs
path_hist = [handles.file_path, '\Track_No_', ...
    num2str(handles.tr_id, '%03d') ,'_', num2str(handles.dp, ...
    '%04d'), 'dp_Speed_hist.tif'];
savePlotWithinGUI(handles.hist_axes, path_hist)

path_speed_time = [handles.file_path, '\Track_No_', ...
    num2str(handles.tr_id, '%03d') ,'_', num2str(handles.dp, ...
    '%04d') ,'dp_Speed_vs_time.tif'];
% speed_vs_time = getframe(handles.speed_vs_time_axes);
% imwrite(speed_vs_time.cdata, path_speed_time, 'tif'); 
savePlotWithinGUI(handles.speed_vs_time_axes, path_speed_time)

path_track = [handles.file_path, '\Track_No_', ...
    num2str(handles.tr_id, '%03d') ,'_', num2str(handles.dp, ...
    '%04d'),'dp_track.tif'];
tr_movie = getframe(handles.movie_axes);
imwrite(tr_movie.cdata, path_track, 'tif'); 

%assign data for single track;
handles.data(:,1) = handles.tr_id;
handles.data(:,2) = handles.time;
handles.data(:,3:4) = handles.pos;
handles.data(:,5) = handles.speed;

%update variables
handles.total_dp = handles.total_dp + handles.dp;
guidata(hObject, handles);

%store data and data_all
path_data = [handles.file_path, '\Track_No_', num2str(handles.tr_id, ...
    '%03d'), '_data.txt'];
data_save = handles.data;                                                           
save(path_data, 'data_save', '-ascii', '-double', '-tabs');


 
function savePlotWithinGUI(axesObject, path)

%axesObject is the axes object that will be saved (required input)

%create a new figure
newFig = figure('Visible', 'off');
 
%set the units and get the position of the axes object
set(axesObject,'Units', 'pixels');
axes_pos = get(axesObject,'Position');

%axes width and height in units of pixels
axes_width = axes_pos(3);
axes_height = axes_pos(4);

%copies axesObject onto new figure
axesObject2 = copyobj(axesObject,newFig);

 
set(newFig, 'Position', [0, 50, axes_width+100, axes_height+100]); 

%realign the axes object on the new figure
set(axesObject2,'Units','pixels');
set(axesObject2,'Position',[80 60 axes_width axes_height]);
set(axesObject2,'Units','normalized');
%saves the plot
saveas(newFig,path, 'tif') 
 
%closes the figure
close(newFig)



function skip_frames_Callback(hObject, eventdata, handles)
% hObject    handle to skip_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of skip_frames as text
%        str2double(get(hObject,'String')) returns contents of skip_frames as a double


% --- Executes during object creation, after setting all properties.
function skip_frames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to skip_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in prev_tr.
function prev_tr_Callback(hObject, eventdata, handles)
% hObject    handle to prev_tr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.tr_id > 1)
    handles.tr_id = handles.tr_id - 1;
end
set(handles.curr_tr, 'String', int2str(handles.tr_id));
guidata(hObject, handles);
evaluate_track(handles);



function curr_tr_Callback(hObject, eventdata, handles)
% hObject    handle to curr_tr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of curr_tr as text
%        str2double(get(hObject,'String')) returns contents of curr_tr as a double

new_tr_number = round(str2double(get(hObject, 'String')));

if ~isnan(new_tr_number)
    if (new_tr_number < 1 || new_tr_number > handles.num_tracks)
        set(hObject, 'String', '1');
        return;
    end
    handles.tr_id = new_tr_number;
    set(handles.curr_tr, 'String', int2str(handles.tr_id));
    guidata(hObject, handles);
    evaluate_track(handles);
end



function zoom_rg_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_rg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zoom_rg as text
%        str2double(get(hObject,'String')) returns contents of zoom_rg as a double

handles.img_zoom = str2double(get(handles.zoom_rg, 'String'));
pos_pxl = handles.pos/handles.magni;              
guidata(hObject, handles);
update_movie(pos_pxl, handles, 1);

% --- Executes during object creation, after setting all properties.
function zoom_rg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zoom_rg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_all.
function save_all_Callback(hObject, eventdata, handles)
% hObject    handle to save_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%assign data_all
%erase null elements for saving
[FileName, PathName] = uigetfile('data.txt', 'Choose first data-file');
if (isscalar(FileName) && isscalar(PathName))
    return;
end
cd(PathName)
data_files = dir('*data.txt');
num_files = size(data_files, 1);

data_all = zeros(1000000,5);
for i=1:num_files
    data = load(data_files(i).name);
    null_ind = find(data_all(:,1)==0,1);
    ind = null_ind:null_ind+size(data,1)-1;
    data_all(ind,:) = data;
end

data_all(find(data_all(:,1)==0,1):end,:)=[];

path_data_all = [handles.file_path, '\data_all.txt'];
save(path_data_all, 'data_all', '-ascii', '-double', '-tabs');

x_hist = 0:1/20:5;
n = hist(data_all(:,5)', x_hist);
n = n./sum(n);
%plot speed histogram
set(handles.hist, 'XData', x_hist, 'YData', n); 

%plot speed versus time
set(handles.speed_time, 'LineStyle', 'none', 'Marker', '+')
set(handles.speed_time, 'XData', data_all(:,2), 'YData', data_all(:,5));
set(handles.time_ind, 'Visible', 'off');
    
%fitting a gaussian to speed histogram
gaussian = 1;
if (gaussian == 1)
    f = fittype('gauss1');
elseif (gaussian == 2)
    f = fittype('gauss2');
end
res_fit  = fit(x_hist(2:end)',n(2:end)',f);
confi_speed = confint(res_fit);

if (gaussian == 1)
     err_speed = confi_speed(2,2) - confi_speed(1,2);
elseif (gaussian == 2)
     err_speed_1 = confi_speed(2,2) - confi_speed(1,2);
     err_speed_2 = confi_speed(2,5) - confi_speed(1,5);
end

hold(handles.hist_axes, 'on')
plot(handles.hist_axes, x_hist, res_fit(x_hist), '-k'),
hold(handles.hist_axes, 'off')
set(handles.evaluation, 'CurrentAxes', handles.hist_axes)
y_lim = get(handles.hist_axes, 'YLim');
y_max = y_lim(end);
if (gaussian == 1)
    text(.1,y_max, ...
    ['Single Gaussian Fit: \newline\Rightarrow v =(',...
        num2str(res_fit.b1,'%1.3f'),'\pm', num2str(err_speed, ...
        '%1.3f)µm/s'), '\newline Mean velocity =', ...
        num2str(mean(data_all(:,5)), '%1.3fµm/s')],...
        'VerticalAlignment', 'top'),
elseif (gaussian == 2)
text(.1,y_max, 'VerticalAlignment', 'top',...
['Double Gaussian Fit\newline\Rightarrow v_1=(',...
    num2str(res_fit.b1,'%1.3f'),'\pm', num2str(err_speed_1, ...
    '%1.3f)µm/s'),'\newline\Rightarrow v_2=(',...
    num2str(res_fit.b2,'%1.3f'),'\pm', num2str(err_speed_2, ...
    '%1.3f)µm/s')]),
end

path_hist = [handles.file_path, '\Speed_all_histogram.tif'];
savePlotWithinGUI(handles.hist_axes, path_hist)
path_speed_time = [handles.file_path, '\Speed_all_vs_time.tif'];
savePlotWithinGUI(handles.speed_vs_time_axes, path_speed_time)



function mov_1st_frame_out_Callback(hObject, eventdata, handles)
% hObject    handle to mov_1st_frame_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mov_1st_frame_out as text
%        str2double(get(hObject,'String')) returns contents of mov_1st_frame_out as a double

mov_1st_frame = round(str2double(get(hObject,'String')));
set(hObject, 'String', int2str(mov_1st_frame));
if mov_1st_frame < 1
    handles.mov_1st_frame = 1;
elseif mov_1st_frame > handles.dp
    handles.mov_1st_frame = handles.dp;
end
pos_pxl = handles.pos/handles.magni;
guidata(handles.evaluation, handles);
update_movie(pos_pxl, handles, mov_1st_frame)

set(hObject, 'String', int2str(handles.mov_1st_frame));
set(handles.frame_slider, 'Value', handles.mov_1st_frame);


% --- Executes during object creation, after setting all properties.
function mov_1st_frame_out_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mov_1st_frame_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function frame_slider_Callback(hObject, eventdata, handles)
% hObject    handle to frame_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

frame = round(get(hObject,'Value'));

pos_pxl = handles.pos/handles.magni;
handles.mov_1st_frame = frame;
guidata(handles.evaluation, handles);
update_movie(pos_pxl, handles, frame)


set(handles.mov_1st_frame_out, 'String', int2str(frame));


% --- Executes during object creation, after setting all properties.
function frame_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in plot_msd.
function plot_msd_Callback(hObject, eventdata, handles)
% hObject    handle to plot_msd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

show_msd = get(handles.plot_msd, 'Value');
if show_msd==1
    set(handles.msd_plot_fig, 'Visible', 'on')
else
    set(handles.msd_plot_fig, 'Visible', 'off');
end

% --- Executes on button press in plot_acf.
function plot_acf_Callback(hObject, eventdata, handles)
% hObject    handle to plot_acf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

show_acf = get(handles.plot_acf, 'Value');
if show_acf==1
    set(handles.acf_plot_fig, 'Visible', 'on')
else
    set(handles.acf_plot_fig, 'Visible', 'off');
end


function delta_t_in_Callback(hObject, eventdata, handles)
% hObject    handle to delta_t_in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delta_t_in as text
%        str2double(get(hObject,'String')) returns contents of delta_t_in as a double

handles.delta_t = str2double(get(hObject,'String'));
handles.dist = round(handles.delta_t*handles.fps);
if handles.dist<1
    handles.dist = 1;
elseif handles.dist>handles.fr_id(end)
    handles.dist = handles.fr_id(end-10);
end
handles.delta_t = handles.dist/handles.fps;
set(hObject, 'String', num2str(handles.delta_t, '%1.1f'));
guidata(handles.evaluation, handles);
evaluate_track(handles)

% --- Executes during object creation, after setting all properties.
function delta_t_in_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delta_t_in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in auto.
function auto_Callback(hObject, eventdata, handles)
% hObject    handle to auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject, 'Value')==1
    for i=1:handles.num_tracks
        handles.tr_id = i;
        set(handles.curr_tr, 'String', int2str(handles.tr_id));
        guidata(hObject, handles);
        evaluate_track(handles);
    end
end
