%% Program just to determine MSD of tracks from Fijji tracks Excel spreadsheet %%
% Program used to determine MSD regimes for colony eversion project,
% Author=Isabelle Wielert


%Format of plots
set(groot, 'DefaultTextInterpreter', 'LaTeX');
set(groot, 'DefaultAxesTickLabelInterpreter', 'LaTeX');
set(groot, 'DefaultAxesFontName', 'LaTeX');
set(groot, 'DefaultLegendInterpreter', 'LaTeX');
set(0,'DefaultAxesColorOrder',summer(16));
%read in .csv data from tracking of Fiji (here without prior processing with stackreg)
string = 'Folder name\';
data = readtable(strcat(string,strcat('duringP6','.csv')));
input = data;




%% bring data of .csv in right format and read out

input=table2array(input);
input(1:3,:)=[]; % delete first three rows (just text and information)

MSD_to_mean=[];
DATA=[];
for i=1:size(input,2)
    
    
    S = sprintf('%s ', input{:,i});
    form= sscanf(S, '%f');
    DATA=[DATA form];
    
end


global num_fr DT trl_t num_tr;
num_fr = max(DATA(:,8)) + 1; % +1 because time_id starts at 0
DT=5; %time steps (Framerate)
trl_t = 8; % threshold length of tracks, default is 3, after that is
T = (0:1:num_fr-1)*DT; k=1; % bins per mum
warning('off','all');
% input=DATA;

%% write MASTER part 1: x,y,vx,vy,COM -------------------------------------
[M,num_tr] = master(DATA);


%% determine MSD
MSD_tracks={};
Rad_pos=[];
Error={};
figure();
for i=1:size(M,1)
    %benötigte Variablen aus Mastermatrix raus lesen
    Rowy=M(i,:,2);
    Rowx=M(i,:,1);
    Rowy(isnan(Rowy)) = [];
    Rowx(isnan(Rowx)) = [];
    Mvelox=M(i,:,3);
    Mveloy=M(i,:,4);
    Mveloy(isnan( Mveloy)) = [];
    Mvelox(isnan( Mvelox)) = [];
    
    
    plot(Rowx,Rowy)
    hold on
    Com=[];
    
    
    meanmsd=zeros(1,length(Rowx));
    N=length(Rowx);
    rad_posi=zeros(1,length(Rowx));
    C_v=zeros(1,length(Rowx));
    
    
    for j=1:N-1 %VACF berechnen
        for k=j:N-1
            C_vhelp=Mvelox(j)*Mvelox(k)+Mveloy(j)*Mveloy(k);
            o = k-j+1;
            C_v(o)=C_v(o)+C_vhelp;
        end
    end
    
    for j=1:N-1 % MSD berechnen zu verschiedenen t
        for k=j+1:N-1
            dist=abs((Rowx(j)-Rowx(k))^2+(Rowy(j)-Rowy(k))^2);
            h = k-j;
            meanmsd(h)=meanmsd(h) + dist;
        end
    end
    
    for j=1:N-1 %über die Zeit mitteln
        meanmsd(j)=meanmsd(j)/(N-j+1);
        C_v(j)=C_v(j)./(N-j+1);
    end
    
    MSD_tracks{1,i}=meanmsd;
    
    Error{1,i}=C_v; %alle VACF in Error speichern
    
end


hold off








%% VACFS average (if you would like to substract that as error)
% Bring data in right format to process (average) it further
Error(:,find(all(cellfun(@isempty, Error),1))) = []; % delete all empty cells in MSD_tracks!
maxSizeEr = max(cellfun(@numel, Error));    %# Get the maximum vector size
fcner = @(x) [x  nan(1,maxSizeEr-numel(x))];  %# semicolon here
rmaterror = cellfun(fcner,Error,'UniformOutput',false);  %# Pad each cell with NaNs
rmaterror=rmaterror';

all_C_v=cell2mat(rmaterror);

nanmean_C_v=nanmean(all_C_v,1);



figure();
plot(nanmean_C_v)
hold off

%% MSD matricses in right format

MSD_tracks(:,find(all(cellfun(@isempty, MSD_tracks),1))) = []; % delete all empty cells in MSD_tracks!
maxSize = max(cellfun(@numel, MSD_tracks));    %# Get the maximum vector size
fcn = @(x) [x  nan(1,maxSize-numel(x))];  %# semicolon here
rmat = cellfun(fcn, MSD_tracks,'UniformOutput',false);  %# Pad each cell with NaNs
rmat=rmat';


allMSD=cell2mat(rmat); % to plot single
plotMSD=cell2mat(rmat); % to mean

%% MSD mitteln über verschiedene RADII (Abstand von Center of mass)
plotMSD=nanmean(plotMSD(:,1:end));
for i=1:maxSize
    plotMSD=[plotMSD nan];
end
MSD_to_mean=[MSD_to_mean; plotMSD];


%% plot different MSD independent of shell

figure();
linspace=[1:maxSize]*DT;
for i=1:length(rmat)
    loglog(linspace(1:end),allMSD(i,:),'r');
    hold on;
end
loglog(linspace(1:end),nanmean(allMSD,1),'b');
hold on;

plot([1:1:60],0.01*[1:1:60],'b');
plot([1:1:60],0.01*[1:1:60].^2,'b');
grid on
hold off

%% saveMSD

filename = strcat(string,strcat('filename','.xlsx'));

