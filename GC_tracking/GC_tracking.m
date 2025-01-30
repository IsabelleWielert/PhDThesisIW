function varargout = GC_tracking(varargin)
% GC_TRACKING M-file for GC_tracking.fig
%      GC_TRACKING, by itself, creates a new GC_TRACKING or raises the existing
%      singleton*.
%
%      H = GC_TRACKING returns the handle to a new GC_TRACKING or the handle to
%      the existing singleton*.
%
%      GC_TRACKING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GC_TRACKING.M with the given input arguments.
%
%      GC_TRACKING('Property','Value',...) creates a new GC_TRACKING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tracking_program_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GC_tracking_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help GC_tracking

% Last Modified by GUIDE v2.5 22-Oct-2011 19:59:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GC_tracking_OpeningFcn, ...
                   'gui_OutputFcn',  @GC_tracking_OutputFcn, ...
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


% --- Executes just before GC_tracking is made visible.
function GC_tracking_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GC_tracking (see VARARGIN)

% Choose default command line output for GC_tracking
handles.output = hObject;

%set graphical objects for histogram
set(0, 'CurrentFigure', handles.GC_tracking);
set(handles.GC_tracking, 'CurrentAxes', handles.img_axes)
colormap(handles.img_axes, 'gray');
img = zeros(512,512);
handles.img = image(img);
axis(handles.img_axes, 'image');
hold(handles.img_axes, 'on')
handles.particles = plot(handles.img_axes, 0, 0, '+r');
hold(handles.img_axes, 'off')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GC_tracking wait for user response (see UIRESUME)
% uiwait(handles.GC_tracking);


% --- Outputs from this function are returned to the command line.
function varargout = GC_tracking_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --- Executes on selection change in camera_input.
function camera_input_Callback(hObject, eventdata, handles)
% hObject    handle to camera_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns camera_input contents as cell array
%        contents{get(hObject,'Value')} returns selected item from camera_input

% Update handles structure
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function camera_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to camera_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function num_folders_Callback(hObject, eventdata, handles)
% hObject    handle to num_folders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_folders as text
%        str2double(get(hObject,'String')) returns contents of num_folders as a double

% Update handles structure
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function num_folders_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_folders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function diameter_input_Callback(hObject, eventdata, handles)
% hObject    handle to diameter_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of diameter_input as text
%        str2double(get(hObject,'String')) returns contents of diameter_input as a double

% Update handles structure
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function diameter_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diameter_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function tracking_on_Callback(hObject, eventdata, handles)
% hObject    handle to tracking_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get user input from GUI

%run tracking routine
track_gonococci(handles);

 
% Update handles structure
guidata(hObject, handles);



% --- Executes on selection change in tif_format.
function tif_format_Callback(hObject, eventdata, handles)
% hObject    handle to tif_format (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns tif_format contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tif_format


% --- Executes during object creation, after setting all properties.
function tif_format_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tif_format (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function binning_Callback(hObject, eventdata, handles)
% hObject    handle to binning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of binning as text
%        str2double(get(hObject,'String')) returns contents of binning as a double


% --- Executes during object creation, after setting all properties.
function binning_CreateFcn(hObject, eventdata, handles)
% hObject    handle to binning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in obj_app.
function obj_app_Callback(hObject, eventdata, handles)
% hObject    handle to obj_app (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns obj_app contents as cell array
%        contents{get(hObject,'Value')} returns selected item from obj_app


% --- Executes during object creation, after setting all properties.
function obj_app_CreateFcn(hObject, eventdata, handles)
% hObject    handle to obj_app (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in eval_tr.
function eval_tr_Callback(hObject, eventdata, handles)
% hObject    handle to eval_tr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

camera = get(handles.camera_input,'Value');
tif_format = get(handles.tif_format,'Value');
binning = str2double(get(handles.binning,'String'));
fps = str2double(get(handles.fps,'String'));
delta_t = str2double(get(handles.delta_t,'String'));


evaluate_tracks(camera, tif_format, binning, fps, delta_t);


function fps_Callback(hObject, eventdata, handles)
% hObject    handle to fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fps as text
%        str2double(get(hObject,'String')) returns contents of fps as a double


% --- Executes during object creation, after setting all properties.
function fps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function delta_t_Callback(hObject, eventdata, handles)
% hObject    handle to delta_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delta_t as text
%        str2double(get(hObject,'String')) returns contents of delta_t as a double


% --- Executes during object creation, after setting all properties.
function delta_t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delta_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function track_gonococci(handles)


obj_app = get(handles.obj_app,'Value');
camera = get(handles.camera_input,'Value');
diameter = str2double(get(handles.diameter_input,'String'));
num_folders = str2double(get(handles.num_folders,'String'));
binning = str2double(get(handles.binning,'String'));
tif_format = get(handles.tif_format,'Value');
thrh_set = str2double(get(handles.thrh, 'String'));
fps = str2double(get(handles.fps,'String'));

%open dialog for choosing the first tif-image in the desired folder
[FileName,PathName] = uigetfile('*.tif; *.bmp');
%if cancel is pushed, then return.
if (isscalar(FileName) == 1) && (isscalar(PathName) == 1);
    return;
end  

% Ending of the image file (e.g. .bmp, .tif);
path_end = FileName(end-3:end);

% pixel factor for absolute lengthscales
switch camera
    case 1 %pixel scaling in units of micron for 
           %100x oil immersion objective (TIRF setup)
       % magni = binning*0.160;
        magni=binning*0.083;
    case 2 %pixel scaling in units of micron for 
           %100x oil immersion objective (Confocal setup)
        magni = binning*0.083;
    case 3  %pixel scaling in units of micron for 
            %100x oil immersion objective (Zeiss OT setup)
        magni = binning*0.099;
    case 4 %pixel scaling in units of micron for
           %100x oil immersion objective (attic lab, OT setup)
        magni = binning*0.059;
end        

%calc diameter in units of pixel
diameter = round(diameter/magni)+1;
if (diameter < 2) || (diameter > 100)
    warndlg('Diameter of object is out of range');
    return;
elseif (diameter/2 == floor(diameter/2))
    diameter = round(diameter/2)*2 + 1;    
end

%find paths of folders
cd(PathName);
folder_path = cell(num_folders,1);
if num_folders == 1;
    folder_path = cellstr(cd);
else
    cd('..')
    Main_Path = cd;
    Main_Path_dir = dir(Main_Path);
    k = 1;
    for z=1:size(Main_Path_dir,1)
       if Main_Path_dir(z).isdir == 1 && strcmp(Main_Path_dir(z).name, '.')==0 ...
           && strcmp(Main_Path_dir(z).name, '..')==0; 
          folder_path(k) = cellstr([Main_Path, '\', Main_Path_dir(z).name]);
          k = k + 1;
       end
    end
end

%id of first good frame
frame_1st = 1;

%Parameters for subroutine track.m
param.mem = round(str2double(get(handles.t_recovery, 'String'))*fps);
                   %this is thenumber of time steps that a particle can be
                   %'lost' and then recovered again.
param.dim = 2; %dimension of the positionlist, here 2 for x- and y-values
param.good = round(str2double(get(handles.min_t, 'String'))*fps);
%              set this keyword to eliminate all trajectories with
%              fewer than param.good valid positions.  This is useful
%              for eliminating very short, mostly 'lost' trajectories
%              due to blinking 'noise' particles in the data stream.
param.quiet = 0;%set to 1 if you want to eliminate any text
max_displ = round(str2double(get(handles.input_max_displ, 'String'))/magni);

for f=1:num_folders
    folder = [char(folder_path(f)),'\'];
    files = dir([folder, '*.tif']);
    
    %delete possible folders in files to get only tiff files
    num_files = size(files,1);
    for j=1:num_files
        if isdir([folder, files(j).name])
            files(j) = [];
        end
    end
    %get image info
    im_info = imfinfo([folder, files(1).name]);
    im_height = im_info.Height;
    im_width = im_info.Width;
    set(handles.img_axes, 'XLim', [1 im_width], 'YLim', [1 im_height]);
    
    %get file name
    if (tif_format <= 2)
        %check if files contain multipage tif file
        if size(im_info, 1) > 1
            set(handles.data_output, 'String', ['Wrong tif format! ',...
            'File is a multipage format']);
            return;
        end
        str = findstr(files(1).name, 't');
        FileName = files(1).name(1:str(1)-1);
        frames = size(files,1);
    else
        FileName = files(1).name(1:end-4);
        frames = size(im_info,1);
    end
    if (isdir([folder, 'results_', FileName])==0)
        mkdir([folder, 'results_', FileName]);   
    end
    
    %define paths for track data
    tbl_file_mat = [folder, 'results_', FileName,'\track_', FileName, '.mat'];
    tbl_file_txt = [folder, 'results_', FileName,'\track_', FileName, '.txt'];
    
    pos_list = zeros(1000000,3);
    %tic to estimate remaining time for tracking
    tic
    %Analyze Frames
    set(handles.data_output, 'String', 'Tracking runs!');
    drawnow
    for i=1:frames        
      switch tif_format
        case 1 %single 8-bit tif files (Camware software for PCO)
            img = double(imread([folder,files(i).name]));            
        case 2 %single 16-bit tif files (NIS Elements at TIRF microscope)
            img = double(imread([folder,files(i).name])); 
        case 3 %multipaged 8-bit tif files (Camware software for PCO)
            img = double(imread([folder,files(1).name],i)); 
        case 4 %multipaged 16-bit tif files (Camware software for PCO)
            img = double(imread([folder,files(1).name],i)); 
      end
      %do contrast enhancement
      int_min = min(img(:));
      int_max = max(img(:));
      if tif_format == 1 || tif_format == 3 
          img = (img-int_min).*(2^8/(int_max-int_min));
      end
      if tif_format == 2 || tif_format == 4 
          img = (img-int_min).*(2^16/(int_max-int_min));
          %conversion to 8-bit
          img = round(img./(2^8));
      end
      
      img_sv = img;
      %invert the 8bit-image if necessary
      if obj_app == 1
         img = 2^8 - img;
      end
      
      %do contrast enhancement
      int_min = min(img(:));
      int_max = max(img(:));
      img = (img-int_min).*(2^8/(int_max-int_min));
      %do real-space bandpass calc to highlight gonococci areas
      num_dim=size(img,3);
      if num_dim > 1 && i==1
          disp('Lena, yet another RGB image!');
      end
      img=img(:,:,1);
      img = bpass(img, .5, diameter + 2);
      %now calc intensity peaks to get positions with pixel acurracy
      max_bright = max(img(:)); 
      thrh =  thrh_set/100*max_bright;
      pk = pkfnd(img, thrh, diameter + 2);
      %now calc centre of intensity to get positions with subpixel acurracy
      cnt = cntrd(img, pk, diameter + 2);
      num_obj = size(cnt,1); %number of tracked objects 
      
      %check tracking parameters at 1st frame 
      if i==frame_1st && f==1
          if num_obj == 0
            frame_1st = frame_1st + 1;
          else
            set(handles.img, 'CData', img)
            set(handles.particles, 'XData', cnt(:,1), 'YData', cnt(:,2))
            %question dialog for stopping tracking or not.
            drawnow
            button = questdlg('Stop image processing to change tracking parameters?', ...
                'Deficient Tracking?', 'Yes', 'No', 'Yes');
            if (strcmp(button, 'Yes') == 1)
                return;
            end
          end
      end
      
      %at each 10th frame show img and tracked particles and calculate
      % remaining time
      if mod(i,10)==0 && num_obj ~= 0;
        t_remain = round(toc/10*(frames-i));
        t_re_min = floor(t_remain/60);
        t_re_sec = mod(t_remain, 60);
        set(handles.data_output, 'String',['Folder: ', int2str(f), ...
            '; Frames to go: ', int2str(frames-i), '; Threshold: ', ...
            num2str(thrh, '%3.2f'), '; Remaining Time: ', ...
            num2str(t_re_min, '%3.0fmin'), ' ', num2str(t_re_sec, '%3.0fs')]);
        set(handles.img, 'CData', img)
        set(handles.particles, 'XData', cnt(:,1), 'YData', cnt(:,2))
        drawnow
        tic
      end
      
      null_id = find(pos_list(:,1)==0,1);
      if num_obj ~= 0
          ind = null_id:null_id + num_obj - 1;
          pos_list(ind, 1:2) = cnt(:,1:2);
          pos_list(ind, 3) = i;
      else
          pos_list(null_id, 3) = i;
      end
      
      % buffer track data each 100th frame
      if mod(i,100)==0
        set(handles.data_output, 'String','Buffering Data:');
        pos_list_buffer = pos_list;
        pos_list_buffer(find(pos_list_buffer(:,1)==0,1):end,:)=[];  
        tr = track(pos_list_buffer, max_displ, handles.data_output, param);
        tr(:,1:2) = magni*tr(:,1:2);
        save(tbl_file_mat, 'tr', '-mat');
        save(tbl_file_txt, 'tr', '-ascii', '-double' ,'-tabs'); 
        drawnow
        tic
        clear pos_list_buffer tr
      end
      
      if get(handles.tracking_off, 'Value') == 1
          if tif_format <= 2
              write_multipage_tif(handles.data_output, folder, files, ...
                                           path_end, tif_format, i, frames)
          end
          set(handles.tracking_off, 'Value', 0)
          return;
      end
      
%       %write multipage tif, if files are single tifs and delete single tifs 
%       if tif_format<=2
%           imwrite(uint8(img_sv), [folder, FileName, path_end], 'tif', ...
%                 'WriteMode', 'append', 'Compression', 'none')
%           delete([folder,files(i).name]);
%       end
          
    end

    pos_list(find(pos_list(:,1)==0,1):end,:)=[];

    tr = track(pos_list, max_displ, handles.data_output, param);
    tr(:,1:2) = magni*tr(:,1:2);
    save(tbl_file_mat, 'tr', '-mat');
    save(tbl_file_txt, 'tr', '-ascii', '-double' ,'-tabs');  
    clear pos_list tr
end



function write_multipage_tif(h_output, folder, files, path_end, tif_format, i, frames)

str = findstr(files(1).name, '_');
FileName = files(1).name(1:str(end)-1);
for j=i:frames
    switch tif_format
        case 1 %single 8-bit tif files (Camware software for PCO)
            % string for new files
            img = double(imread([folder,files(j).name]));            
        case 2 %single 16-bit tif files (NIS Elements at TIRF microscope)
            img = double(imread([folder,files(j).name])); 
            %conversion to 8-bit
            img = round(img./(2^8));
    end
    imwrite(uint8(img), [folder, FileName, path_end], 'tif', ...
                    'WriteMode', 'append', 'Compression', 'none')
    delete([folder,files(j).name]);
    if mod(j,10)==0
        set(h_output, 'String',['Complete generation of ',...
            'multipage tif of remaining files. Frames to go: ',...
            int2str(frames-j)]);
        drawnow
    end
end

% --- Executes on button press in tracking_off.
function tracking_off_Callback(hObject, eventdata, handles)
% hObject    handle to tracking_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function thrh_Callback(hObject, eventdata, handles)
% hObject    handle to thrh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thrh as text
%        str2double(get(hObject,'String')) returns contents of thrh as a double

thrh = round(str2double(get(hObject,'String')));
set(hObject, 'String', int2str(thrh));
if thrh > 100
    thrh = 100;
    set(hObject, 'String', int2str(thrh));
elseif thrh < 1
    thrh = 1;
    set(hObject, 'String', int2str(thrh));
end
set(handles.thrh_slider, 'Value', thrh);

% --- Executes during object creation, after setting all properties.
function thrh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thrh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in test_thrh.
function test_thrh_Callback(hObject, eventdata, handles)
% hObject    handle to test_thrh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (get(hObject, 'Value') == 1);
    obj_app = get(handles.obj_app,'Value');
    camera = get(handles.camera_input,'Value');
    diameter = str2double(get(handles.diameter_input,'String'));
    binning = str2double(get(handles.binning,'String'));
    tif_format = get(handles.tif_format,'Value');

    %open dialog for choosing the first tif-image in the desired folder
    [FileName,PathName] = uigetfile('*.tif; *.bmp');
    %if cancel is pushed, then return.
    if (isscalar(FileName) == 1) && (isscalar(PathName) == 1);
        return;
    end  

    % pixel factor for absolute lengthscales
    switch camera
        case 1 %pixel scaling in units of micron for 
               %100x oil immersion objective (TIRF setup)
            %magni = binning*0.160;
            magni = binning*0.083;
        case 2 %pixel scaling in units of micron for 
               %100x oil immersion objective (Confocal setup)
            magni = binning*0.0833;
        case 3  %pixel scaling in units of micron for 
                %100x oil immersion objective (Zeiss OT setup)
            magni = binning*0.099;
        case 4 %pixel scaling in units of micron for
               %100x oil immersion objective (attic lab, OT setup)
            magni = binning*0.059;
    end        

    %calc diameter in units of pixel
    diameter = round(diameter/magni)+1;
    if (diameter < 2) || (diameter > 100)
        warndlg('Diameter of object is out of range');
        return;
    elseif (diameter/2 == floor(diameter/2))
        diameter = round(diameter/2)*2 + 1;    
    end
    cd(PathName);
    file_path = [PathName, FileName];

    %get image info
    im_info = imfinfo(file_path);
    im_height = im_info.Height;
    im_width = im_info.Width;
    set(handles.img_axes, 'XLim', [1 im_width], 'YLim', [1 im_height]);
    drawnow

    %get file name
    if (tif_format <= 2)
        %check if files contain multipage tif file
        if size(im_info, 1) > 1
            set(handles.data_output, 'String', ['Wrong tif format! ',...
            'File is a multipage format']);
            return;
        end
    end

    switch tif_format
        case 1 %single 8-bit tif files (Camware software for PCO)
            img = double(imread(file_path));            
        case 2 %single 16-bit tif files (NIS Elements at TIRF microscope)
            img = double(imread(file_path)); 
        case 3 %multipaged 8-bit tif files (Camware software for PCO)
            img = double(imread(file_path,1)); 
        case 4 %multipaged 16-bit tif files (Camware software for PCO)
            img = double(imread(file_path,1)); 
    end
    %do contrast enhancement
    int_min = min(img(:));
    int_max = max(img(:));
    if tif_format == 1 || tif_format == 3 
        img = (img-int_min).*(2^8/(int_max-int_min));
    end
    if tif_format == 2 || tif_format == 4 
        img = (img-int_min).*(2^16/(int_max-int_min));
        %conversion to 8-bit
        img = round(img./(2^8));
    end

    %invert the 8bit-image if necessary
    if obj_app == 1
     img = 2^8 - img;
    end
    
    %do contrast enhancement
    int_min = min(min(img));
    int_max = max(max(img));
    img_ct = (img-int_min).*(2^8/(int_max-int_min));
    %do real-space bandpass calc to highlight gonococci areas
    img_filter = bpass(img_ct, .5, diameter + 2);
    %now calc intensity peaks to get positions with pixel acurracy
    max_bright = max(max(img_filter)); 
    set(handles.img, 'CData', img_filter)

    while (get(hObject, 'Value') == 1)
        thrh = str2double(get(handles.thrh, 'String'));
        thrh =  thrh/100*max_bright;
        pk = pkfnd(img_filter, thrh, diameter + 2);
        %now calc centre of intensity to get positions with subpixel acurracy
        cnt = cntrd(img_filter, pk, diameter + 2); 
        if isempty(cnt)
            h = warndlg(['There were no positions inputted into calculations.',...
                          '. check your Intensity threshold!'], '!!Warning!!');
            uiwait(h);   
            set(handles.thrh, 'String', '50');
            thrh_Callback(handles.thrh, eventdata, handles)
            cnt = [-1,-1];
        end
        set(handles.particles, 'XData', cnt(:,1), 'YData', cnt(:,2))
        drawnow
        pause(0.5);
    end    
end

% --- Executes on slider movement.
function thrh_slider_Callback(hObject, eventdata, handles)
% hObject    handle to thrh_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

thrh = get(hObject,'Value');
set(handles.thrh, 'String', int2str(thrh));

% --- Executes during object creation, after setting all properties.
function thrh_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thrh_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function min_t_Callback(hObject, eventdata, handles)
% hObject    handle to min_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_t as text
%        str2double(get(hObject,'String')) returns contents of min_t as a double


% --- Executes during object creation, after setting all properties.
function min_t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t_recovery_Callback(hObject, eventdata, handles)
% hObject    handle to t_recovery (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_recovery as text
%        str2double(get(hObject,'String')) returns contents of t_recovery as a double


% --- Executes during object creation, after setting all properties.
function t_recovery_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_recovery (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function input_max_displ_Callback(hObject, eventdata, handles)
% hObject    handle to input_max_displ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_max_displ as text
%        str2double(get(hObject,'String')) returns contents of input_max_displ as a double


% --- Executes during object creation, after setting all properties.
function input_max_displ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_max_displ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
