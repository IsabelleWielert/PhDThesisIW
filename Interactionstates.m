function varargout = Interactionstates(varargin)
% INTERACTIONSTATES MATLAB code for Interactionstates.fig
%      INTERACTIONSTATES, by itself, creates a new INTERACTIONSTATES or raises the existing
%      singleton*.
%
%      H = INTERACTIONSTATES returns the handle to a new INTERACTIONSTATES or the handle to
%      the existing singleton*.
%
%      INTERACTIONSTATES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTERACTIONSTATES.M with the given input arguments.
%
%      INTERACTIONSTATES('Property','Value',...) creates a new INTERACTIONSTATES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Interactionstates_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Interactionstates_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Interactionstates

% Last Modified by GUIDE v2.5 30-Jan-2020 15:33:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Interactionstates_OpeningFcn, ...
                   'gui_OutputFcn',  @Interactionstates_OutputFcn, ...
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


% --- Executes just before Interactionstates is made visible.
function Interactionstates_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Interactionstates (see VARARGIN)

% Choose default command line output for Interactionstates
handles.output = hObject;
handles.datatrackoben=[];
handles.datatrackunten=[];
handles.fullpath=[];
handles.fullpath2=[];
handles.velocities=[];
handles.force=[];
handles.Fullcaltime=[];
handles.start=[];
handles.end=[];
handles.times=[];

handles.x=[];

handles.current_data1=[];
handles.Cal_time=[];
handles.probability2=[];
handles.probability1=[];
handles.probability3=[];
handles.probability4=[];
handles.probability5=[];
handles.data1=[];
handles.data2=[];
handles.data3=[];
handles.data4=[];
handles.data5=[];

handles.dataforce1=[];
handles.dataforce2=[];
handles.k1= 0.013; %pN/nm
handles.k2= 0.012; %pN/nm

handles.framerate=50;
% handles.Masstab=106.83;%RobertsFalle
 handles.Masstab=88.652;%MeineFalle -> Überprüft
% handles.Masstab=58;
handles.counter1 = 0;
handles.counter2 = 0;
handles.counter3 = 0;
handles.counter4 = 0;
handles.counter5 = 0;
handles.counterforce=0;
handles.counterforce2=0;
handles.d=5; %Zeit für downsampling der GESCHWINDIGKEIT
axes(handles.Trackoben)
xlabel(handles.Trackoben,'Time [s]')
ylabel(handles.Trackoben,'Displacement [nm]')

axes(handles.Trackunten)
xlabel(handles.Trackunten,'Time [s]')
ylabel(handles.Trackunten,'Displacement[nm]')


axes(handles.ForceVelocity)

datacursormode on

yyaxis(handles.ForceVelocity,'left');
set(handles.ForceVelocity, 'YTickMode', 'auto')
handles.plot_force = plot(0,0,'-k',0,0,'-k');
YLabel = get(handles.ForceVelocity, 'YLabel');
set(YLabel, 'String', 'force');

yyaxis(handles.ForceVelocity,'right');
handles.plot_velocity = plot(0,0,'-r',0,0,'-r','LineWidth',0.5);
set(handles.ForceVelocity, 'YTickMode', 'auto')
YLabel = get(handles.ForceVelocity, 'YLabel');
set(YLabel, 'String', 'velocity');
ylim([-5000 2000]);

XLabel = get(handles.ForceVelocity, 'XLabel');
set(XLabel, 'String', 'time [s]');
Title = get(handles.ForceVelocity, 'Title');
set(Title, 'String', 'Force and Velocity', 'FontWeight', 'bold');
grid(handles.ForceVelocity,'on');
set(handles.ForceVelocity,'XMinorGrid','on','XMinorTick','on');


%set paraneter für tables
set(handles.table_elo,'Data',[]);
set(handles.table_retract,'Data',[]);
set(handles.table_paus,'Data',[]);
set(handles.table_bundl,'Data',[]);
set(handles.table_no,'Data',[]);
set(handles.table_rupt,'Data',[]);
set(handles.table_revers,'Data',[]);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Interactionstates wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Interactionstates_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Path1.
function Path1_Callback(hObject, eventdata, handles)
% hObject    handle to Path1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile('C:\Users\Isabelle\Desktop\ToAnalysenewpilELISA\RichtigerTail17\ToanalyseData\Day6\*.txt');
if isequal(file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(path,file)]);
end
foldername=strsplit(file,'.');
handles.namesforsave=foldername{1};
handles.fullpath=fullfile(path,file);
% 
%    fid2=fopen(handles.fullpath);
%    handles.data_oben=fread(fid2,[2,inf],'double','b');
%    fclose(fid2);
 handles.data_oben=textread(handles.fullpath);
% handles.data_oben=flip(handles.data_oben);
handles.data_oben=handles.data_oben';
handles.data_veloben=handles.data_oben;

handles.data_oben=handles.data_oben.*handles.Masstab;
handles.x=[1:length(handles.data_oben(1,1:end))]./handles.framerate;
%   handles.x=[1:length(handles.data_oben(1,1:end))]./handles.framerate;%Robertsdaten
handles.times=handles.x;
% handles.x=[1:length(handles.data_oben(2,1:end))]./handles.framerate;
% handles.times=handles.x;
handles.DeltaT=handles.d/handles.framerate;
   
 
axes(handles.Trackoben)
 plot(handles.x,handles.data_oben(1,1:end));
% plot(handles.x,handles.data_oben(1,1:end));%Robetsdaten

set(handles.Pathtext1,'String',handles.fullpath)
guidata(hObject, handles);


% --- Executes on button press in Path2.
function Path2_Callback(hObject, eventdata, handles)
% hObject    handle to Path2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path2] = uigetfile('C:\Users\Isabelle\Desktop\ToAnalysenewpilELISA\RichtigerTail17\ToanalyseData\Day6\*.txt');
if isequal(file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(path2,file)]);
end

handles.fullpath2=fullfile(path2,file);

  % fid2=fopen(handles.fullpath2);
   %handles.data_unten=fread(fid2,[2,inf],'double','b');
   handles.data_unten=textread(handles.fullpath2);
   
   %fclose(fid2);

handles.data_unten=handles.data_unten';
handles.data_velunten=handles.data_unten;
handles.data_unten=handles.data_unten.*handles.Masstab;

 handles.x=[1:length(handles.data_unten(1,1:end))]./handles.framerate;
 % handles.x=[1:length(handles.data_unten(1,1:end))]./handles.framerate;%Robertsdaten
handles.times=handles.x;
handles.DeltaT=handles.d/handles.framerate;
axes(handles.Trackunten)
 plot(handles.x,handles.data_unten(1,1:end));
%plot(handles.x,handles.data_unten(1,1:end));%Robertsdaten

set(handles.Pathtext2,'String',handles.fullpath2)
guidata(hObject, handles);


% --- Executes on button press in Correcttrack1.
function Correcttrack1_Callback(hObject, eventdata, handles)
% hObject    handle to Correcttrack1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.Trackoben)

    h=datacursormode(gcf);
    set(h,'Displaystyle','Window','Enable','on');

while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 1)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x1=dat.Position(1)*handles.framerate;
    y1=dat.Position(2);
    
   
    break
    end
    end
    
    while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 2)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x2=dat.Position(1)*handles.framerate;
    y2=dat.Position(2);
    break
    end
    end
    while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 2)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x3=dat.Position(1)*handles.framerate;
    y3=dat.Position(2);
    break
    end
    end
    while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 2)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(  gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x4=dat.Position(1)*handles.framerate;
    y4=dat.Position(2);
    break
    end
    end
    %Steigung 
    Xfit2=[1:length(handles.data_oben(1,1:end))];
    Xfit2=Xfit2./handles.framerate;
    fit_func = fittype('poly1'); 
    fitdata2 = fit([Xfit2(x1:x2) Xfit2(x3:x4)]',[handles.data_oben(1,x1:x2) handles.data_oben(1,x3:x4)]'./handles.Masstab ,fit_func);
    m1=coeffvalues(fitdata2);
    
    
    handles.baselinesubstraction2=Xfit2.*m1(1)./handles.Masstab+mean(handles.data_oben(1,x1:x2));
    handles.data_oben1 = handles.data_oben(1,1:end)-handles.baselinesubstraction2;
    handles.data_obenplot=handles.data_oben1;
    handles.data_oben=handles.data_obenplot;
    handles.datatrackoben=handles.data_oben;
    plot(handles.x,handles.data_obenplot);

guidata(hObject, handles);

% --- Executes on button press in Correcttrack2.
function Correcttrack2_Callback(hObject, eventdata, handles)
% hObject    handle to Correcttrack2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.Trackunten)
   h=datacursormode(gcf);
    set(h,'Displaystyle','Window','Enable','on');

while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 1)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x1=dat.Position(1)*handles.framerate;
    y1=dat.Position(2);
    
   
    break
    end
    end
    
    while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 2)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x2=dat.Position(1)*handles.framerate;
    y2=dat.Position(2);
    break
    end
    end
    while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 2)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x3=dat.Position(1)*handles.framerate;
    y3=dat.Position(2);
    break
    end
    end
    while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 2)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(  gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x4=dat.Position(1)*handles.framerate;
    y4=dat.Position(2);
    break
    end
    end
    %Steigung 
    Xfit=[1:length(handles.data_unten(1,1:end))];
    Xfit=Xfit./handles.framerate;
    fit_func = fittype('poly1'); 
    fitdata = fit([Xfit(x1:x2) Xfit(x3:x4)]',[handles.data_unten(1,x1:x2) handles.data_unten(1,x3:x4)]'./handles.Masstab ,fit_func);
    m=coeffvalues(fitdata);
    
    
    handles.baselinesubstraction=Xfit.*m(1)./handles.Masstab+mean(handles.data_unten(1,x1:x2));
    handles.data_unten1 = handles.data_unten(1,1:end)-handles.baselinesubstraction;
    handles.data_untenplot=handles.data_unten1;
    handles.data_unten=handles.data_untenplot;
    handles.datatrackunten=handles.data_unten;
    a=length(handles.data_untenplot);
    b=length(handles.data_obenplot);
      plot(handles.x,handles.data_untenplot);
    if a<b
        c=b-a;
        handles.data_oben(end-c+1:end)=[];
        length( handles.data_oben)
        length( handles.data_unten)
        handles.x=[1:length(handles.data_unten(2:end))]./handles.framerate;
    else
         c=a-b;
        handles.data_unten(end-c+1:end)=[];
        length( handles.data_oben)
        length( handles.data_unten)
        handles.x=[1:length(handles.data_oben(2:end))]./handles.framerate;
    end
        
handles.times=handles.x;


guidata(hObject, handles);

% --- Executes on button press in Duration.
function Duration_Callback(hObject, eventdata, handles)
% hObject    handle to Duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Cal_time=length(handles.data_oben)./handles.framerate;
set(handles.showduration,'String',num2str(handles.Cal_time))
handles.Fullcaltime=handles.Cal_time;

handles.data_oben=handles.data_oben;
handles.data_unten=handles.data_unten;


% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.velocity=zeros(1,ceil(length(handles.data_oben(1:handles.d:end))));
j=0;
for i=1:handles.d:ceil(length(handles.data_unten(1:end)))-2*handles.d
    j=j+1;
    handles.velocity(j)=((mean(handles.data_unten(i+handles.d:i+2*handles.d))-mean(handles.data_oben(i+handles.d:i+2*handles.d)))-(mean(handles.data_unten(i:i+handles.d))-mean(handles.data_oben(i:i+handles.d))))/handles.DeltaT;
 if handles.velocity(j)<-2500
     handles.velocity(j)=0
 end
end




handles.force=(handles.data_oben(1:end).*handles.k1 +(handles.data_unten(1:end)).*handles.k2.*(-1));
handles.current_data1=handles.force;
%handles.current_data2=handles.velocity;
%handles.velocities=handles.velocity;
 handles.current_data2=handles.velocity;
 handles.velocities=smooth(handles.velocity);
%handles.TABLEtimeforce=[handles.times, handles.force]



axes(handles.ForceVelocity)

set(handles.plot_force, 'XData',handles.x, 'YData', handles.current_data1(1:length(handles.x)));
% set(handles.plot_velocity, 'XData',handles.x(1:handles.d:end)+handles.d/handles.framerate, 'YData', handles.current_data2);
set(handles.plot_velocity, 'XData',handles.x(1:handles.d:end)+handles.d/handles.framerate, 'YData', handles.current_data2);

% set(handles.plot_velocity, 'XData',handles.x(1:handles.d:end)+handles.d/handles.framerate, 'YData', handles.current_data2);
   h=datacursormode(gcf);
    set(h,'Displaystyle','Window','Enable','on'); 
    xlim([0 handles.Cal_time]);

guidata(hObject, handles);

% --- Executes on button press in Cut.
function Cut_Callback(hObject, eventdata, handles)
axes(handles.ForceVelocity)
   h=datacursormode(gcf);
    set(h,'Displaystyle','Window','Enable','on');

% hObject    handle to Cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 1)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x1=dat.Position(1)*handles.framerate;

    
   
    break
    end
    end
    
    while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 2)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x2=dat.Position(1)*handles.framerate;

    break
    end
    end
    

    %Steigung 
    handles.start=x1;
    handles.end=x2;
    
handles.x=handles.x(x1:x2)-x1/handles.framerate;
handles.Cal_time=length(handles.x)./handles.framerate
handles.current_data1=handles.current_data1(x1:x2);
handles.current_data2=handles.current_data2(ceil(x1/handles.d):ceil(x2/handles.d));
handles.velocity=handles.current_data2;
langeint=length(handles.x(1:handles.d:end));
handles.current_data2=handles.current_data2(1:langeint)
set(handles.plot_force, 'XData',handles.x, 'YData', handles.current_data1);
set(handles.plot_velocity, 'XData',handles.x(1:handles.d:end)+handles.d/handles.framerate, 'YData', handles.current_data2);
% set(handles.plot_velocity, 'XData',handles.x(1:handles.d:end)+handles.d/handles.framerate, 'YData', handles.current_data2);

handles.xvel=handles.x(1:handles.d:end)+handles.d/handles.framerate;

set(handles.showduration,'String',num2str(handles.Cal_time))
xlim([0 handles.Cal_time]);

guidata(hObject, handles);

% --- Executes on button press in analyse.
function analyse_Callback(hObject, eventdata, handles)
% hObject    handle to analyse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.ForceVelocity)
yyaxis(handles.ForceVelocity,'right');
yline(0);
yline(110);
yline(-110);

yyaxis(handles.ForceVelocity,'left')
ylim([-5 20]);
xlim([0 11]);

guidata(hObject, handles);

% --- Executes on button press in prob_retr.
function prob_retr_Callback(hObject, eventdata, handles)
% hObject    handle to prob_retr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datatable=get(handles.table_retract,'Data');
duration=0;
for i = 1: length(datatable(:,1))
duration= duration+abs(datatable{i,1}-datatable{i,2})/handles.framerate;
end
handles.probability2=duration/handles.Cal_time
set(handles.textprob_retr,'String',num2str(handles.probability2))


guidata(hObject, handles);

% --- Executes on button press in prob_elo.
function prob_elo_Callback(hObject, eventdata, handles)
% hObject    handle to prob_elo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datatable=get(handles.table_elo,'Data');
duration=0;
for i = 1: length(datatable(:,1))
duration= duration+abs(datatable{i,1}-datatable{i,2})/handles.framerate;
end
handles.probability1=duration/handles.Cal_time
set(handles.textprob_elo,'String',num2str(handles.probability1))


guidata(hObject, handles);

% --- Executes on button press in prob_paus.
function prob_paus_Callback(hObject, eventdata, handles)
% hObject    handle to prob_paus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datatable=get(handles.table_paus,'Data');
duration=0;
for i = 1: length(datatable(:,1))
duration= duration+abs(datatable{i,1}-datatable{i,2})/handles.framerate;
end
handles.probability3=duration/handles.Cal_time
set(handles.textprob_paus,'String',num2str(handles.probability3))


guidata(hObject, handles);

% --- Executes on button press in prob_bundl.
function prob_bundl_Callback(hObject, eventdata, handles)
% hObject    handle to prob_bundl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datatable=get(handles.table_bundl,'Data');
duration=0;
for i = 1: length(datatable(:,1))
duration= duration+abs(datatable{i,1}-datatable{i,2})/handles.framerate;
end
handles.probability4=duration/handles.Cal_time
set(handles.textprob_bundl,'String',num2str(handles.probability4))


guidata(hObject, handles);

% --- Executes on button press in prob_no.
function prob_no_Callback(hObject, eventdata, handles)
% hObject    handle to prob_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datatable=get(handles.table_no,'Data');
duration=0;
for i = 1: length(datatable(:,1))
duration= duration+abs(datatable{i,1}-datatable{i,2})/handles.framerate;
end
handles.probability5=duration/handles.Cal_time
set(handles.textprob_no,'String',num2str(handles.probability5))


guidata(hObject, handles);


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

framerate=handles.framerate;
stiffnes1=handles.k1;
stiffnes2=handles.k2;
path1=handles.fullpath;
path2=handles.fullpath2;
trackobenfull=handles.datatrackoben;
trackuntenfull=handles.datatrackunten;
velocityfull=handles.velocities;
forcefull=handles.force;
durationtime=handles.Fullcaltime;
cuttingstart=handles.start;
cuttingend=handles.end;
timesfull=handles.times;

cuttedtimes=handles.x;

cuttedforce=handles.current_data1;
durationcutted=handles.Cal_time;
cuttedprobrup=handles.probability2;
cuttedprobelo=handles.probability1;
cuttedprobpaus=handles.probability3;
cuttedprobbundl=handles.probability4;
cuttedprobnoev=handles.probability5;
Rupturestate=handles.data1;
Elongationstate=handles.data2;
Pausingstate=handles.data3;
Bundlingstate=handles.data4;
NoEventsstate=handles.data5;

Ruptureforces=handles.dataforce1;

Reversalforces=handles.dataforce2;


folderrawAuswertung = ['C:\Users\Isabelle\Desktop\ToAnalysenewpilELISA\RichtigerTail17\Analyse10pr\day6_analysis\',handles.namesforsave,];
% folderrawAuswertung = ['C:\Users\Isabelle\Desktop\Lisa_Day1and2Auswertung\Toanalyse_Lisas_day1and2\G4_analysed\day1_analysis\',handles.Dateiname.String,];handles.namesforsave
if exist(folderrawAuswertung, 'dir') == 0
      mkdir(folderrawAuswertung)
end



savefile = [folderrawAuswertung,'\Auswertung.mat'];

 
save(savefile,'framerate','stiffnes2','stiffnes2','path1','path2','trackobenfull','trackuntenfull','timesfull','velocityfull','forcefull','durationtime','cuttingstart','cuttingend','cuttedtimes','cuttedforce','durationcutted','cuttedprobrup','cuttedprobelo','cuttedprobpaus','cuttedprobbundl','cuttedprobnoev','Rupturestate','Elongationstate','Pausingstate','Bundlingstate','NoEventsstate','Ruptureforces','Reversalforces');



guidata(hObject, handles);
    

function Dateiname_Callback(hObject, eventdata, handles)
% hObject    handle to Dateiname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dateiname as text
%        str2double(get(hObject,'String')) returns contents of Dateiname as a double


% --- Executes during object creation, after setting all properties.
function Dateiname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dateiname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in retraction.
function retraction_Callback(hObject, eventdata, handles)
% hObject    handle to retraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.ForceVelocity)
h=datacursormode(gcf);
set(h,'Displaystyle','Window','Enable','on');
 
while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 1)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x1=dat.Position(1)*handles.framerate;
    y1=dat.Position(2)*handles.framerate;
    
   
    break
    end
    end
    
    while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 2)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x2=dat.Position(1)*handles.framerate;
    y2=dat.Position(2)*handles.framerate;
    break
    end
    end
    
    
    meanvel=mean(handles.velocity(floor((x1)/handles.d):floor(x2/(handles.d))));


handles.counter1 = handles.counter1 + 1;
handles.data1 = get(handles.table_retract,'Data');

 if isempty(handles.data1)
     handles.data1 = {x1 x2 meanvel};
  else
% %     % append the new row

handles.data1 = [handles.data1; {x1 x2 meanvel}];

 end
  set(handles.table_retract,'Data',handles.data1);
  handles.data1 = get(handles.table_retract,'Data');
%plot(handles.x,handles.current_data1,'k*')
hold on 


for i=1:handles.counter1
plot(handles.x(handles.data1{i,1}:handles.data1{i,2}),handles.current_data1(handles.data1{i,1}:handles.data1{i,2}),'-c')
x1=handles.data1{end,1}/handles.framerate;
hold on
end

pan on
guidata(hObject, handles);


% --- Executes on button press in pausing.
function pausing_Callback(hObject, eventdata, handles)
% hObject    handle to pausing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.ForceVelocity)
h=datacursormode(gcf);
set(h,'Displaystyle','Window','Enable','on');

while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 1)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x1=dat.Position(1)*handles.framerate;
    y1=dat.Position(2)*handles.framerate;
    
   
    break
    end
    end
    
    while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 2)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x2=dat.Position(1)*handles.framerate;
    y2=dat.Position(2)*handles.framerate;
    break
    end
    end

handles.counter3 = handles.counter3 + 1;
handles.data3 = get(handles.table_paus,'Data');

 if isempty(handles.data3)
     handles.data3 = {x1 x2};
  else
% %     % append the new row

handles.data3 = [handles.data3; {x1 x2}];

 end
  set(handles.table_paus,'Data',handles.data3);
  handles.data3 = get(handles.table_paus,'Data');
hold on 
for i=1:handles.counter3
plot(handles.x(handles.data3{i,1}:handles.data3{i,2}),handles.current_data1(handles.data3{i,1}:handles.data3{i,2}),'-b')
x1=handles.data3{end,1}/handles.framerate;
hold on
end
pan on
guidata(hObject, handles);

% --- Executes on button press in noevent.
function noevent_Callback(hObject, eventdata, handles)
% hObject    handle to noevent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.ForceVelocity)
h=datacursormode(gcf);
set(h,'Displaystyle','Window','Enable','on');
while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 1)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x1=dat.Position(1)*handles.framerate;
    y1=dat.Position(2)*handles.framerate;
    
   
    break
    end
    end
    
    while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 2)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x2=dat.Position(1)*handles.framerate;
    y2=dat.Position(2)*handles.framerate;
    break
    end
    end
    


handles.counter5 = handles.counter5 + 1;
handles.data5 = get(handles.table_no,'Data');

 if isempty(handles.data5)
     handles.data5 = {x1 x2};
  else
% %     % append the new row

handles.data5 = [handles.data5; {x1 x2}];

 end
  set(handles.table_no,'Data',handles.data5);
  handles.data5 = get(handles.table_no,'Data');

hold on 
for i=1:handles.counter5
plot(handles.x(handles.data5{i,1}:handles.data5{i,2}),handles.current_data1(handles.data5{i,1}:handles.data5{i,2}),'-y')
x1=handles.data5{end,1}/handles.framerate;

hold on
end
pan on
guidata(hObject, handles);



% --- Executes on button press in bundling.
function bundling_Callback(hObject, eventdata, handles)
% hObject    handle to bundling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.ForceVelocity)
h=datacursormode(gcf);
set(h,'Displaystyle','Window','Enable','on');
while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 1)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x1=dat.Position(1)*handles.framerate;
    y1=dat.Position(2)*handles.framerate;
    
   
    break
    end
    end
    
    while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 2)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x2=dat.Position(1)*handles.framerate;
    y2=dat.Position(2)*handles.framerate;
    break
    end
    end
    


handles.counter4 = handles.counter4 + 1;
handles.data4 = get(handles.table_bundl,'Data');

 if isempty(handles.data4)
     handles.data4 = {x1 x2};
  else
% %     % append the new row

handles.data4 = [handles.data4; {x1 x2}];

 end
  set(handles.table_bundl,'Data',handles.data4);
  handles.data4 = get(handles.table_bundl,'Data');
hold on 
for i=1:handles.counter4
plot(handles.x(handles.data4{i,1}:handles.data4{i,2}),handles.current_data1(handles.data4{i,1}:handles.data4{i,2}),'-m')
x1=handles.data4{end,1}/handles.framerate;
hold on
end
pan on
guidata(hObject, handles);

% --- Executes on button press in elongation.
function elongation_Callback(hObject, eventdata, handles)
% hObject    handle to elongation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.ForceVelocity)
h=datacursormode(gcf);
set(h,'Displaystyle','Window','Enable','on');
while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 1)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x1=dat.Position(1)*handles.framerate;
    y1=dat.Position(2)*handles.framerate;
    
   
    break
    end
    end
    
    while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 2)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x2=dat.Position(1)*handles.framerate;
    y2=dat.Position(2)*handles.framerate;
    break
    end
    end
meanvel=mean(handles.velocity(ceil((x1+3)/handles.d):handles.d:ceil((x2-5)/handles.d)));

handles.counter2 = handles.counter2 + 1;
handles.data2 = get(handles.table_elo,'Data');

 if isempty(handles.data2)
     handles.data2 = {x1 x2 meanvel};
  else
% %     % append the new row

handles.data2 = [handles.data2; {x1 x2 meanvel}];

 end
  set(handles.table_elo,'Data',handles.data2);
  handles.data2 = get(handles.table_elo,'Data');
hold on 
for i=1:handles.counter2
plot(handles.x(handles.data2{i,1}:handles.data2{i,2}),handles.current_data1(handles.data2{i,1}:handles.data2{i,2}),'-r')
x1=handles.data2{end,1}/handles.framerate;
hold on
end
pan on
guidata(hObject, handles);




% --- Executes on button press in rupture.
function rupture_Callback(hObject, eventdata, handles)
% hObject    handle to rupture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.ForceVelocity)
   h=datacursormode(gcf);
    set(h,'Displaystyle','Window','Enable','on');
while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 1)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x1=dat.Position(1);
    y1=dat.Position(2);
    
   
    break
    end
end
handles.dataforce1 = get(handles.table_rupt,'Data');

 if isempty(handles.dataforce1)
     handles.dataforce1 = {x1 y1};
  else
% %     % append the new row

handles.dataforce1 = [handles.dataforce1; {x1 y1}];

 end
  set(handles.table_rupt,'Data',handles.dataforce1);
  handles.dataforce1 = get(handles.table_rupt,'Data');
pan on
guidata(hObject, handles);
    

% --- Executes on button press in reversal.
function reversal_Callback(hObject, eventdata, handles)
% hObject    handle to reversal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.ForceVelocity)
h=datacursormode(gcf);
set(h,'Displaystyle','Window','Enable','on');
while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 1)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x1=dat.Position(1);
    y1=dat.Position(2);
    
   
    break
    end
end
handles.dataforce2 = get(handles.table_revers,'Data');

 if isempty(handles.dataforce2)
     handles.dataforce2 = {x1 y1};
  else
% %     % append the new row

handles.dataforce2 = [handles.dataforce2; {x1 y1}];

 end
  set(handles.table_revers,'Data',handles.dataforce2);
  handles.dataforce2 = get(handles.table_revers,'Data');

pan on
guidata(hObject, handles);
