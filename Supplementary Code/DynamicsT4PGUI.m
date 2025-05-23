function varargout = Dynamicsbearb180222(varargin)
% DYNAMICSBEARB180222 MATLAB code for Dynamicsbearb180222.fig
%      DYNAMICSBEARB180222, by itself, creates a new DYNAMICSBEARB180222 or raises the existing
%      singleton*.
%
%      H = DYNAMICSBEARB180222 returns the handle to a new DYNAMICSBEARB180222 or the handle to
%      the existing singleton*.
%
%      DYNAMICSBEARB180222('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DYNAMICSBEARB180222.M with the given input arguments.
%
%      DYNAMICSBEARB180222('Property','Value',...) creates a new DYNAMICSBEARB180222 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Dynamicsbearb180222_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Dynamicsbearb180222_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Dynamicsbearb180222

% Last Modified by GUIDE v2.5 18-Feb-2022 16:50:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Dynamicsbearb180222_OpeningFcn, ...
                   'gui_OutputFcn',  @Dynamicsbearb180222_OutputFcn, ...
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


% --- Executes just before Dynamicsbearb180222 is made visible.
function Dynamicsbearb180222_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Dynamicsbearb180222 (see VARARGIN)

% Choose default command line output for Dynamicsbearb180222
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
handles.current_data1=[];
handles.Cal_time=[];
handles.data1=[];
handles.data2=[];
handles.data3=[];
handles.data4=[];
handles.data5=[];
handles.piluslengths=[];
handles.dataforce1=[];
handles.dataforce2=[];
ginfo=load('info.mat')
handles.arr2=ginfo.arr2;
handles.framerate=ginfo.framerate;
handles.savepath=ginfo.allsavepath
 handles.Masstab= 0.0792354;%Orca-> Überprüft

handles.counter1 = 0;
handles.counter2 = 0;
handles.counter3 = 0;
handles.counter4 = 0;
handles.counter5 = 0;
handles.counterforce=0;
handles.counterforce2=0;
g=load('track.mat')
 handles.data=[g.XDATA g.YDATA];
 
 

axes(handles.ForceVelocity)
datacursormode on

set(handles.ForceVelocity, 'YTickMode', 'auto')
handles.plot_force = plot(handles.data(:,1),handles.data(:,2),'k');
hold on 
YLabel = get(handles.ForceVelocity, 'YLabel');
set(YLabel, 'String', 'Length [Px]');


XLabel = get(handles.ForceVelocity, 'XLabel');
set(XLabel, 'String', 'Frames');
Title = get(handles.ForceVelocity, 'Title');
set(Title, 'String', 'Length of Pilus', 'FontWeight', 'bold');
grid(handles.ForceVelocity,'on');
set(handles.ForceVelocity,'XMinorGrid','on','XMinorTick','on');


%set paraneter für tables
set(handles.table_elo,'Data',[]);
set(handles.table_retract,'Data',[]);
set(handles.table_paus,'Data',[]);
set(handles.LengthofPilus,'Data',[]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Dynamicsbearb180222 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Dynamicsbearb180222_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Duration.
function Duration_Callback(hObject, eventdata, handles)
% hObject    handle to Duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles.x=[1:length(handles.data)]./handles.framerate
 handles.lengthi=max(handles.data);

handles.Cal_time=handles.lengthi(1)./handles.framerate;
set(handles.showduration,'String',num2str(handles.Cal_time))
handles.Fullcaltime=handles.Cal_time;



axes(handles.ForceVelocity)
h=datacursormode(gcf);
set(h,'Displaystyle','Window','Enable','on'); 
hold off
    
guidata(hObject, handles);



% --- Executes on button press in Cut.
function Cut_Callback(hObject, eventdata, handles)
% hObject    handle to Cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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
    x1=dat.Position(1);
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
    x2=dat.Position(1);
    y2=dat.Position(2);

    break
    end
    end
    


    
ToBeFoundOut=[x1 y1];
ToBeFoundOut2=[x2 y2];


L = find(all(bsxfun(@eq,handles.data',ToBeFoundOut')));

L2 = find(all(bsxfun(@eq,handles.data',ToBeFoundOut2')));

if L2<L
handles.x=handles.data(L2:L,1);
handles.y=handles.data(L2:L,2);
else
handles.x=handles.data(L:L2,1);
handles.y=handles.data(L:L2,2);
end
    
    
    
% handles.datalen=handles.data(x1:x2,1);
handles.datalen=handles.x;
handles.Cal_time=length(handles.datalen)./handles.framerate;
% handles.x=handles.data(length(handles.data)-x2-y1:length(handles.data)-x1-y1,1);
% handles.y=handles.data(length(handles.data)-x2-y1:length(handles.data)-x1-y1,2);

handles.Arrayana=[handles.x handles.y];
set(handles.plot_force, 'XData',handles.x, 'YData', handles.y);

set(handles.showduration,'String',num2str(handles.Cal_time))


guidata(hObject, handles);





% --- Executes on button press in Elongation.
function Elongation_Callback(hObject, eventdata, handles)
% hObject    handle to Elongation (see GCBO)
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
    
    while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 2)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x2=dat.Position(1);
    y2=dat.Position(2);
    break
    end
    end
    
% meanvel=((y2-y1)/(x2-x1))*handles.framerate*handles.Masstab;
ToBeFoundOut=[x1 y1];
ToBeFoundOut2=[x2 y2];


L = find(all(bsxfun(@eq,handles.Arrayana',ToBeFoundOut')));

L2 = find(all(bsxfun(@eq,handles.Arrayana',ToBeFoundOut2')));



[xData,yData]=prepareCurveData([L:L2],handles.y(L:L2));

ft=fittype('poly1');
opts=fitoptions('Method','LinearLeastSquares');
[fitresult,gof]=fit(xData, yData, ft, opts);
% figure(); 
% h=plot(fitresult,xData,yData)
c=coeffvalues(fitresult);

meanvel=c(1)*handles.framerate*handles.Masstab;
hold on 
plot(handles.x(L:L2),c(1).*[L:L2]+c(2))
hold on 
% 
% velocitiees=[];
% for s = [L:L2-1]
%     veloci=((handles.y(s+1)-handles.y(s))/1)*handles.framerate*handles.Masstab;
% velocitiees=[velocitiees;veloci];
% end
% meanvel=mean(veloci);

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


ToBeFoundOut=[x1 y1];
ToBeFoundOut2=[x2 y2];


L = find(all(bsxfun(@eq,handles.Arrayana',ToBeFoundOut')));

L2 = find(all(bsxfun(@eq,handles.Arrayana',ToBeFoundOut2')));

for i=1:handles.counter2
    if L2<L
plot(handles.x(L2:L),handles.y(L2:L),'-r')
    else
plot(handles.x(L:L2),handles.y(L:L2),'-r')
    end
hold on
end

pan on
guidata(hObject, handles);


% --- Executes on button press in Retraction.
function Retraction_Callback(hObject, eventdata, handles)
% hObject    handle to Retraction (see GCBO)
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
    
    while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 2)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x2=dat.Position(1);
    y2=dat.Position(2);
    break
    end
    end
    
% meanvel=((y2-y1)/(x2-x1))*handles.framerate*handles.Masstab;

ToBeFoundOut=[x1 y1];
ToBeFoundOut2=[x2 y2];


L = find(all(bsxfun(@eq,handles.Arrayana',ToBeFoundOut')));

L2 = find(all(bsxfun(@eq,handles.Arrayana',ToBeFoundOut2')));



[xData,yData]=prepareCurveData([L:L2],handles.y(L:L2));

ft=fittype('poly1');
opts=fitoptions('Method','LinearLeastSquares');
[fitresult,gof]=fit(xData, yData, ft, opts);
% figure(); 
% h=plot(fitresult,xData,yData)
c=coeffvalues(fitresult);

meanvel=c(1)*handles.framerate*handles.Masstab;
hold on 
plot(handles.x(L:L2),c(1).*[L:L2]+c(2))
hold on 

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
hold on 
ToBeFoundOut=[x1, y1];
ToBeFoundOut2=[x2, y2];


ToBeFoundOut=[x1 y1];
ToBeFoundOut2=[x2 y2];


L = find(all(bsxfun(@eq,handles.Arrayana',ToBeFoundOut')));

L2 = find(all(bsxfun(@eq,handles.Arrayana',ToBeFoundOut2')));

for i=1:handles.counter1
    if L2<L
plot(handles.x(L2:L),handles.y(L2:L),'-c')
    else
plot(handles.x(L:L2),handles.y(L:L2),'-c')
    end
hold on
end

pan on
guidata(hObject, handles);

% --- Executes on button press in Pausing.
function Pausing_Callback(hObject, eventdata, handles)
% hObject    handle to Pausing (see GCBO)
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
    
    while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 2)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    x2=dat.Position(1);
    y2=dat.Position(2);
    break
    end
    end
    
% meanvel=((y2-y1)/(x2-x1))*handles.framerate*handles.Masstab;


ToBeFoundOut=[x1 y1];
ToBeFoundOut2=[x2 y2];


L = find(all(bsxfun(@eq,handles.Arrayana',ToBeFoundOut')));

L2 = find(all(bsxfun(@eq,handles.Arrayana',ToBeFoundOut2')));



[xData,yData]=prepareCurveData([L:L2],handles.y(L:L2));

ft=fittype('poly1');
opts=fitoptions('Method','LinearLeastSquares');
[fitresult,gof]=fit(xData, yData, ft, opts);
% figure(); 
% h=plot(fitresult,xData,yData)
c=coeffvalues(fitresult);

meanvel=c(1)*handles.framerate*handles.Masstab;
hold on 
plot(handles.x(L:L2),c(1).*[L:L2]+c(2))
hold on 


handles.counter3 = handles.counter3 + 1;
handles.data2 = get(handles.table_paus,'Data');

 if isempty(handles.data3)
     handles.data3 = {x1 x2 meanvel};
  else
% %     % append the new row

handles.data3 = [handles.data3; {x1 x2 meanvel}];

 end
  set(handles.table_paus,'Data',handles.data3);
  handles.data3 = get(handles.table_paus,'Data');
hold on 


ToBeFoundOut=[x1 y1];
ToBeFoundOut2=[x2 y2];


L = find(all(bsxfun(@eq,handles.Arrayana',ToBeFoundOut')));

L2 = find(all(bsxfun(@eq,handles.Arrayana',ToBeFoundOut2')));

for i=1:handles.counter2
    if L2<L
plot(handles.x(L2:L),handles.y(L2:L),'-b')
    else
plot(handles.x(L:L2),handles.y(L:L2),'-b')
    end
hold on
end

pan on
guidata(hObject, handles);



% --- Executes on button press in length.
function length_Callback(hObject, eventdata, handles)
% hObject    handle to length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to Pausing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% prompt9 = {'How many pili lengths do you want to measure ? '};
% dlgtitle9 = 'Input';
% dims9 = [1 100];
% definput9 = {'3'};
% answer9 = inputdlg(prompt9,dlgtitle9,dims9,definput9)  
% Input9=cellfun(@(x)str2double(x), answer9);
% Pilinumber=Input9(1);

axes(handles.ForceVelocity)
h=datacursormode(gcf);
set(h,'Displaystyle','Window','Enable','on');
%  for j = 1:Pilinumber
while 1 %%%%%%%%%%%%%%%%%%%%% LOOP to measure velocity function (POINT 1)
         control=waitforbuttonpress; % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    if control==1
    cursorobj = datacursormode(gcf);
    cursorobj.SnapToDataVertex = 'on';
    dat=getCursorInfo(cursorobj);
    y1=dat.Position(2)*handles.Masstab;
    
   
    break
    end
end
    

handles.piluslengths=[handles.piluslengths; y1]

  
handles.counter4 = handles.counter4 + 1;
handles.data4 = get(handles.LengthofPilus,'Data');

 if isempty(handles.data4)
     handles.data4 = {y1};
  else
% %     % append the new row

handles.data4 = [handles.data4; {y1}];

 end
  set(handles.LengthofPilus,'Data',handles.data4);
  handles.data4 = get(handles.LengthofPilus,'Data');
hold on 
%  end

pan on
guidata(hObject, handles);


% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


pixel=handles.y;

cuttedtimes=handles.x;

durationcutted=handles.Cal_time;

Retractionstate=handles.data1;
Elongationstate=handles.data2;
Pausingstate=handles.data3;
Lengths=handles.data4;


 savepath=cell2mat(handles.savepath);
%  arr2=cell2mat(handles.arr2);
  savefilename = [savepath,'\',handles.arr2];
 savefile = fullfile([savefilename,'\',handles.edit2.String,'.mat'])

% save(filename, 'Pili' );

save(savefile,'Retractionstate','Elongationstate','Pausingstate','Lengths','pixel','cuttedtimes','durationcutted');





guidata(hObject, handles);



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit2.
function edit2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
