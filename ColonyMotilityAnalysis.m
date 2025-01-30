% This program is important to analyse the radial motility of colonies (used for the data:)
clear all;
clc;
close all;

%Track videos of colony (Methods described in Wielert et al. doi:
%https://doi.org/10.1101/2023.07.06.548055) via Fiji (Trackmate PlugIn)
% Import .csv file of tracking:
input = readtable('your file.csv');

% Import one frame of the analysed video (.tif) to determine the center of mass
allpath = Tiff('Your file.tif');
savepath='Your Path';
%% Center of mass
xyscale = 0.0792354;
% Read in image, which was also used for tracking in fiji to
% determine center of mass (COM)
img=read(allpath);
allmin = min(img(:));
allminhelp=allmin;
allmax = max(img(:));
allmaxhelp=allmax;
img(1,1,1) = allmin;
img=uint8(fix(((1500.0+1)*single(img-allmin)-1)/single(allmax-allmin)));
%% Determine center of mass
%Image needs to be pre-processed interms f bianrizining and
%contrast enhancement to be able to read out center of mass
%automatically

Moviebinarize=imbinarize(imadjust(imgaussfilt(img(:,:,1),5),[0.9 1]),'global');
figure();
imshow(Moviebinarize)
hold off
CoM=regionprops(Moviebinarize,'centroid'); %regionprops determines center of mass of binarized images
s=CoM; %COM
comglobal=s(1).Centroid.*xyscale; % COM in right scale






%% read out data of .csv file:

input=table2array(input);
input(1:3,:)=[]; % delte first three rows (just information)

MSD_to_mean=[];
DATA=[];
for i=1:size(input,2)
    S = sprintf('%s ', input{:,i});
    form= sscanf(S, '%f');
    DATA=[DATA form];
    
end


% these paramertes needs to be adapted for your data
global num_fr DT trl_t num_tr;
num_fr = max(DATA(:,8)) + 1; % +1 because time_id starts at 0
DT=0.1; %Frameraze
trl_t = 30; % threshold length of tracks, default is 300
T = (0:1:num_fr-1)*DT; k=1; % bins per mum
warning('off','all');
input=DATA;

%% write MASTER part 1: x,y,vx,vy,COM -------------------------------------
% Bring the data in the right format (the master function was developed by Marc Hennes and is reused here)
[M,num_tr] = master(DATA);


%% Calculate mean squared displacements of tracks:
MSD_tracks={};
Rad_pos=[];
Error={};
for i=1:size(M,1)
    %Varibales needed from master data
    Rowy=M(i,:,2);
    Rowx=M(i,:,1);
    Rowy(isnan(Rowy)) = [];
    Rowx(isnan(Rowx)) = [];
    Mvelox=M(i,:,3);
    Mveloy=M(i,:,4);
    Mveloy(isnan( Mveloy)) = [];
    Mvelox(isnan( Mvelox)) = [];
    
    
    Com=[];
    
    
    meanmsd=zeros(1,length(Rowx));
    N=length(Rowx);
    C_v=zeros(1,length(Rowx));
    
    
    %determine center of mass for each track to determine position
    for p =1:length(Rowx)
        x_center=comglobal(1)-Rowx(p);
        y_center=comglobal(2)-Rowy(p);
        CenterofMass=sqrt(x_center^2+y_center^2); % Radial COM
        Com=[Com;CenterofMass];
    end
    Rad_pos=[Rad_pos;mean(Com)]; % Collect all COM
    binpos=mean(Com); %Mean Position der Tracks
    
    
    
    %determine velocities from master (M(:,:,4&5))
    Velox=zeros(1, length(Rowx)-1);
    Veloy=zeros(1, length(Rowx)-1);
    Radialrowx=comglobal(1)-Rowx;
    Radialrowy=comglobal(2)-Rowy;
    for ti=1:length(Radialrowx)-1
        Velox(ti)=(Radialrowx(ti+1)-Radialrowx(ti))/DT;
        Veloy(ti)=(Radialrowy(ti+1)-Radialrowy(ti))/DT;
    end
    
    
    for j=1:N-1 %determine velocity auto correlation function (VACF)
        for k=j:N-1
            C_vhelp=Velox(j)*Velox(k)+Veloy(j)*Veloy(k);
            o = k-j+1;
            C_v(o)=C_v(o)+C_vhelp;
        end
    end
    
    for j=1:N-1 % determine mean squared displacements
        for k=j+1:N-1
            dist=abs((Rowx(j)-Rowx(k))^2+(Rowy(j)-Rowy(k))^2);
            h = k-j;
            meanmsd(h)=meanmsd(h) + dist;
        end
    end
    for j=1:N-1 %average over time
        meanmsd(j)=meanmsd(j)/(N-j+1);
        C_v(j)=C_v(j)./(N-j+1);
    end
    
    if C_v(2)<0 % Error of tracks are first negative Peakf of VACF (if it is negative, might be a check for your data)
        MSD_tracks{1,i}=meanmsd-2*(DT)^2.*C_v(2); %Fehler abziehen!
    else
        MSD_tracks{1,i}=meanmsd;
    end
    Error{1,i}=C_v; %save in Error
    
end
hold off;

%% Density of the shells (if there too less cells, the analysis is not as conlcusive)
%You can also add a threshhold
Bin=histogram(Rad_pos,'BinWidth',1);
Comparebins=Bin.Values;
Ratio=zeros(1,length(Comparebins));
Lessdense=length(Comparebins);
% for i=1:length(Comparebins)-1
%     Ratio(i)=Comparebins(i+1)/Comparebins(i);
%     Lessdense=i;
%     if Ratio(i)<0.8 && Lessdense >8
%         break
%     end
%
% end


%% VACFS (take the mean to get the overall error if you like)
% Über alle VACFs mitteln
Error(:,find(all(cellfun(@isempty, Error),1))) = []; % delete all empty cells in MSD_tracks!
maxSizeEr = max(cellfun(@numel, Error));    %# Get the maximum vector size
fcner = @(x) [x  nan(1,maxSizeEr-numel(x))];  %# semicolon here
rmaterror = cellfun(fcner,Error,'UniformOutput',false);  %# Pad each cell with NaNs
rmaterror=rmaterror';

all_C_v=cell2mat(rmaterror);

nanmean_C_v=nanmean(all_C_v,1);


% Plot to check how the VACF looks like:
% figure();
% plot(nanmean_C_v)
% hold off

%% MSD Matrices processing (to be able to average them)

MSD_tracks(:,find(all(cellfun(@isempty, MSD_tracks),1))) = []; % delete all empty cells in MSD_tracks!
maxSize = max(cellfun(@numel, MSD_tracks));    %# Get the maximum vector size
fcn = @(x) [x  nan(1,maxSize-numel(x))];  %# semicolon here
rmat = cellfun(fcn, MSD_tracks,'UniformOutput',false);  %# Pad each cell with NaNs
rmat=rmat';


allMSD=cell2mat(rmat); % to plot single
plotMSD=cell2mat(rmat); % to mean

%% MSD averaging depnendent on the distance to the center of mass
plotMSD=nanmean(plotMSD(:,1:end));
for i=1:maxSize
    plotMSD=[plotMSD nan];
end
MSD_to_mean=[MSD_to_mean; plotMSD];


%% plot different MSD independent of shell
% f1=figure();
% figure(f1);
linspace=[1:maxSize]*DT;
% for i=1:length(rmat)
%     loglog(linspace(1:end),allMSD(i,:),'r');
%     hold on;
% end
% % loglog(linspace(1:end),MSD_to_mean(1:maxSize),'b');
% hold off

%% plot different MSD dependent of shell
num_bins_radii=round(max(Rad_pos));
for j = 2:1:num_bins_radii+1
    for i=1:length(Rad_pos)
        if Rad_pos(i)<=j && Rad_pos(i)>j-1
            Cellmsd{j,i}=allMSD(i,:);
        end
    end
end


% Calculate MSD and mean MSD in every shell
for i=1:1:num_bins_radii
    figure();
    Rowi=Cellmsd(i,:);
    Rowi=Rowi';
    Rowmat=cell2mat(Rowi);
    for o=1:size(Rowmat,1)
        jupi=loglog(linspace,Rowmat(o,:),'Color',[0.5 0.5 0.5]);
        jupi.HandleVisibility='off';
        hold on;
    end
    
    MEANMSDGESAMT{i}=nanmean(Rowmat);
    loglog(linspace,MEANMSDGESAMT{i},'r');

    hold off
end


%% Plot Mean MSD in different shells with fit
figure();
for i=1:1:num_bins_radii-1
    if ~isempty(MEANMSDGESAMT{i})
        loglog(linspace,MEANMSDGESAMT{i},'Linewidth',2);
        plot(linspace,MEANMSDGESAMT{i},'Linewidth',2);
    end
    hold on
end
title('Meansquaredisplacement for within colony motility')
x0=200;
y0=100;
width=1350;
height=800;
set(gcf,'position',[x0,y0,width,height]);
legendCell = strcat(string(num2cell([1:num_bins_radii])),'µm');
leg=legend(legendCell,'Location','southeast', 'AutoUpdate', 'off');
title(leg,'Distance to center of colony');
xlabel('lag time [s]')
ylabel('MSD [µm^2]')
hold on
%% Fit for correlation time


hasNaN = cellfun(@nnz,cellfun(@isnan, MEANMSDGESAMT, 'Unif',0), 'Unif',0);    % Cells With ‘NaN’ Values
idx = find([hasNaN{:}]);
for i=idx
    MEANMSDGESAMT{i}={};
end
for i=1:num_bins_radii-1
    if ~isempty(MEANMSDGESAMT{i})
        loglinspace=log(linspace);
        logMSDmean=log(MEANMSDGESAMT{i});
        [xData, yData] = prepareCurveData(loglinspace(1:10),logMSDmean(1:10)); % fit for the first 10 time points, adjust to your data and frame rate
        
        % Set up fittype and options.
        ft = fittype( 'a*x + b', 'independent', 'x', 'dependent', 'y' ); % logarithmic data to perform linear fit
        
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts );
        coefs=coeffvalues(fitresult);
        
        hold on
        yAchse{i}=exp(coefs(2))/4;
        %D=MSD(1s)/4
        %log(4Dt ^alpha)=log(4D)+alpha log(t) -> log(4D)=yAchse
        
        
        
        alpha = 0.95;
        ci{i} = confint(fitresult, alpha) ;
        uncert{i}=coefs(2)-ci{i}(1,2); % Error, this is why VACF is not needed
        ErryAchse{i}=yAchse{i}*uncert{i};
        Exponent{i}=coefs(1);
        Radii{i}=i;
    end
end
hold off

XA=cell2mat(Radii);
YA=cell2mat(yAchse);

ErrorYA=cell2mat(ErryAchse);
save(savepath,'XA','YA','ErrorYA')
%% Figure to plot Diffsuion consant dependent on radial position of shell
figure();
errorbar(XA,YA,ErrorYA);
xlabel('Radii from center [µm]');
ylabel('Diffusion constant [µm^2/s]');