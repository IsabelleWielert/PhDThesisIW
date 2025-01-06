% Script: Determine radial motility along a colony from fluroescent images
% name of script from developmental stage:
% RadialMSD_Isaversionwithoutshells



clear all;
clc;
close all;

%%Formatierung Plots

set(groot, 'DefaultTextInterpreter', 'LaTeX');
set(groot, 'DefaultAxesTickLabelInterpreter', 'LaTeX');
set(groot, 'DefaultAxesFontName', 'LaTeX');
set(groot, 'DefaultLegendInterpreter', 'LaTeX');
set(0,'DefaultAxesColorOrder',summer(7));
%%
%for Marcs Data: change DT and maxSize

%Einlesen der Daten .csv von Fiji, HERE YOU NEED TO ENTER YOUR FILE NAMES

string = 'FOLDER\';
name='NAMEofTrackFilefromTrackmate' %This is a csv file from ImageJ (Fiji) from the Trackmate Plugin, in which the image was processed (contrast) and the single cells were tracked with Trackmate (PlugIn)
nameimage='NameofImage' %This is a .tif file format of one image of the analysed colony to determine the center of mass and the radii

data = readtable(strcat(string,strcat(name,'.csv')));

input = data;

savepath='SaveFolder\';
pathi=[savepath '\Nameoffinalfig.fig']
pathdata=[savepath '\NameofAnalysis.mat']



%Form of data

input=table2array(input);
input(1:3,:)=[]; % ersten drei Textzeilen der Tabelle löschen

MSD_to_mean=[];
DATA=[];
for i=1:size(input,2)
    
    
    S = sprintf('%s ', input{:,i});
    form= sscanf(S, '%f');
    DATA=[DATA form];
    
end


global num_fr DT trl_t num_tr;
num_fr = max(DATA(:,8)) + 1; % +1 because time_id starts at 0
DT=4; %Zeitschritte (Framerate)
trl_t = 3; % threshold length of tracks, default is 3, after that is
T = (0:1:num_fr-1)*DT; k=1; % bins per mum
warning('off','all');
input=DATA;

%Masterdata erstellen (master function bring data in a form to work with (this function was developed by Marc Hennes))
[M,num_tr] = master(DATA);
%% Boundary bestimmen:
xyscale = 0.0792354;

t = Tiff(strcat(string,strcat(nameimage,'.tif')),'r');
imgData = read(t);
% figure();
% imshow(imgData);
% hold off
allmin = min(imgData(:));
allminhelp=allmin;
allmax = max(imgData(:));
allmaxhelp=allmax;
imgData(1,1,1) = allmin;
img=uint8(fix(((2000.0+1)*single(imgData-allmin)-1)/single(allmax-allmin)));
% imshow(img);
bw = imbinarize(imgaussfilt(img,4));
% figure(); imshow(bw); hold off
BW2 = bwareaopen(bw, 1000);
% figure(); imshow(BW2); hold off

se = strel("disk",2);

bw3 = imclose(BW2,se);
% figure(); imshow(bw3); hold off

bw4 = imfill(bw3,"holes");
% figure(); imshow(bw4); hold off

figure();
[B,L] = bwboundaries(bw4,"noholes");
imshow(label2rgb(L,@jet,[.5 .5 .5]))

hold on
for k = 1:length(B)
    boundary = B{k};
    
    plot(boundary(:,2),boundary(:,1),'w','LineWidth',2)
    
    boundaryscaled=boundary.*xyscale;
end
title("Objects with Boundaries in White")


stats = regionprops(L,"Circularity","Centroid", "MajorAxisLength","MinorAxisLength");
threshold = 0;
for k = 1:length(B)
    
    % Obtain the circularity corresponding to label "k"
    circ_value = stats(k).Circularity;
    
    
    % Display the results
    circ_string = sprintf('%2.2f',circ_value);
    
    % Mark objects above the threshold with a black circle
    if circ_value > threshold
        centroid = stats(k).Centroid;
        diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
        
        radiis = (diameters/2)*xyscale;
        plot(centroid(1),centroid(2),'ko');
        hold on
        viscircles(centroid,radiis/xyscale)
    end
    
    text(boundary(1,2)-35,boundary(1,1)+13,circ_string,'Color','y',...
        'FontSize',14,'FontWeight','bold')
    
end
title("Centroids of Circular Objects and Circularity Values")


% If you figure out a problem cause colonies are too large, then you can
% solve it now, just a dialogue if you want to set the center of mass on
% your own:
answer = questdlg('Are you happy with center and radius?', ...
    'Skipping',...
    'yes','no','no');
% Handle response
switch answer
    case 'yes'
        dessert = 1;
    case 'no'
        
        dessert = 2;
end

if dessert == 1
    
    comglobal=[centroid(2) centroid(1)].*xyscale
    
elseif dessert == 2
    
    
    % Set center of mass manually
    f1=figure('Name','Where is the center of the colony?');
    for i=1:size(M,1)
        %benötigte Variablen aus Mastermatrix raus lesen
        Rowy=M(i,:,2);
        Rowx=M(i,:,1);
        plot(Rowx,Rowy,'r') % you get a figure with all the tracks
        
        hold on
    end
    
    coord=ginput(1);hold on
    coord=[coord(2) coord(1)]%.*xyscale
    comglobal=coord
end

% Here a control figure will be plotted

f1=figure('Name','Where is the center of the colony?');
for i=1:size(M,1)
    %benötigte Variablen aus Mastermatrix raus lesen
    Rowy=M(i,:,2);
    Rowx=M(i,:,1);
    plot(Rowx,Rowy,'r')
    hold on
    plot(comglobal(2),comglobal(1),'gx')
    hold on
end


%
%% MSD bestimmen:
MSD_tracks={};
Rad_pos=[];
Error={};
for i=1:size(M,1)
    %benötigte Variablen aus Mastermatrix raus lesen
    Rowy=M(i,:,2);
    Rowx=M(i,:,1);
    Rowy(isnan(Rowy)) = [];
    Rowx(isnan(Rowx)) = [];
    
    plot(Rowx,Rowy,'r')
    
    Mvelox=M(i,:,3);
    Mveloy=M(i,:,4);
    Mveloy(isnan( Mveloy)) = [];
    Mvelox(isnan( Mvelox)) = [];
    
    
    Com=[];
    
    
    meanmsd=zeros(1,length(Rowx));
    N=length(Rowx);
    rad_posi=zeros(1,length(Rowx));
    C_v=zeros(1,length(Rowx));
    
    
    %Center of mass of each track bestimmen um Position
    for p =1:length(Rowx)
        x_center=comglobal(2)-Rowx(p);
        y_center=comglobal(1)-Rowy(p);
        CenterofMass=sqrt(x_center^2+y_center^2); % Radialer Center of Mass
        Com=[Com;CenterofMass];
    end
    Rad_pos=[Rad_pos;mean(Com)]; % Sammel alle CoM
    binpos=mean(Com); %Mean Position der Tracks
    
    
    
    %velocities berechnen by handy from (M(:,:,4&5))
    Velox=zeros(1, length(Rowx)-1);
    Veloy=zeros(1, length(Rowx)-1);
    %     Radialrowx=comglobal(2).-Rowx;
    %     Radialrowy=Rowy;
    Radialrowx=comglobal(2)-Rowx;
    Radialrowy=comglobal(1)-Rowy;
    for ti=1:length(Radialrowx)-1
        Velox(ti)=(Radialrowx(ti+1)-Radialrowx(ti))/DT;
        Veloy(ti)=(Radialrowy(ti+1)-Radialrowy(ti))/DT;
    end
    
    
    for j=1:N-1 %VACF berechnen
        for k=j:N-1
            C_vhelp=Velox(j)*Velox(k)+Veloy(j)*Veloy(k);
            o = k-j+1;
            C_v(o)=C_v(o)+C_vhelp;
        end
    end
    %
    for j=1:N-1 %VACF berechnen
        for k=j:N-1
            C_vhelp=Mvelox(j)*Mvelox(k)+Mveloy(j)*Mveloy(k);
            o = k-j+1;
            C_v(o)=C_v(o)+C_vhelp;
        end
    end
    %
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
    
    %if C_v(2)<0 %Fehler der Tracks ist erster negativer Peak der VACF, nur wenn erster Peak negativ ist
    
    %     MSD_tracks{1,i}=meanmsd-2*(DT)^2.*C_v(2); %Fehler abziehen!
    % else
    MSD_tracks{1,i}=meanmsd;
    %
    Error{1,i}=C_v; %alle VACF in Error speichern
    %
end
hold off;





%% MSD Matrizen vorbereiten
% MSD Matrizen vorbereiten um über alle Tracks mitteln zu können:
maxSize = max(cellfun(@numel, MSD_tracks));
fcn = @(x) [x  nan(1,maxSize-numel(x))];  %# semicolon here
rmat = cellfun(fcn, MSD_tracks,'UniformOutput',false);  %# Pad each cell with NaNs
rmat=rmat';


allMSD=cell2mat(rmat); % to plot single
plotMSD=cell2mat(rmat); % to mean



%% Determine MSD

linspace=[1:maxSize]*DT;


for i=1:length(rmat)
    
    if ~isempty(allMSD(i,:)) && sum(~isnan(allMSD(i,:)) ) >12
        loglinspace=linspace;
        logMSDmean=allMSD(i,:);
        [xData, yData] = prepareCurveData(loglinspace(1:10),logMSDmean(1:10));
        
        % Set up fittype and options.
        ft = fittype( 'a*x + b', 'independent', 'x', 'dependent', 'y' );
        
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        
        % Fit model to data5.
        [fitresult, gof] = fit( xData, yData, ft, opts );
        coefs=coeffvalues(fitresult);
        
        
        yAchse{i}=coefs(1)/4;
        alpha = 0.95;
        ci{i} = confint(fitresult, alpha) ;
        uncert{i}=coefs(1)-ci{i}(1,1);
        ErryAchse{i}=(yAchse{i}/4)*uncert{i};
        Exponent{i}=coefs(2);
        Radii{i}=Rad_pos(i);
    end
end
hold off

XA=cell2mat(Radii);
YA=cell2mat(yAchse);
ErrorYA=cell2mat(ErryAchse);



%% Figure to plot Diffsuion consant dependent on radial position of shell
h=figure();
scatter(XA(1:end),YA(1:end));
ylim([0 0.1])
xlim([1 44])
xlabel('Distance to edge [µm]');
ylabel('Diffusion constant [µm^2/s]');
hold off
%%

[XA_sorted, a_order] = sort(XA,'descend');
XA_sortedfromedge=max(XA_sorted)-XA_sorted;
newYA =YA(a_order);
newYA_error=ErrorYA(a_order)
h3=figure();
plot(XA_sorted(1:end),newYA(1:end),'-r');
ylim([0 0.1])
xlim([0 44])
xlabel('Distance to edge [µm]');
ylabel('Diffusion constant [µm^2/s]');
hold on

M=movmean(newYA,2,"omitnan")
% plot(XA_sorted(1:end),smoothdata(newYA,'movmean',35),'-b');

%% Save the data
save(pathdata,'XA_sortedfromedge','newYA','ErrorYA')
