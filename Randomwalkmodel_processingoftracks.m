%% Program to analyse the tracks of the twichting gonococci and determine velocity and persistence time
% tracking was performed via GC-tracking, which was developed by Lena
% Dewenter
clear;
close all;
scale=0.083;%micron xyscale
DeltaT=1/10; %weil mit 10 Hz gemessen wurde
files=dir('folder od DataAll of Tracking\*.txt');
names={files.name};
PosID=[];
count=1;

for i=1:length(names)
    
    fileID = fopen(['folder od DataAll of Tracking\', names{i}]);
    formatSpec = '%f %f %f %f %f';
    sizeA = [5 ,Inf];
    A = fscanf(fileID,formatSpec,sizeA);
    
    count=count+500;
    datasorth=A';
    
    PosID=[PosID; datasorth(:,3) datasorth(:,4) datasorth(:,1)+count]; % datasorth(:,1)+
    
end

%% Determine the mean squared displacements
MSD_tracks={};

for i=1:max(PosID(:,3))-1
    
    px=PosID(PosID(:,3)==i,1);
    py=PosID(PosID(:,3)==i,2);
    
    meanmsd=zeros(1,length(px)+1);
    N=length(px);
    for j=1:N
        for k=j+1:N
            dist=abs((px(j)-px(k))^2+(py(j)-py(k))^2);
            h = k-j;
            meanmsd(h)=meanmsd(h) + dist;
        end
        
    end
    for j=1:N
        meanmsd(j)=meanmsd(j)/(N-j+1);
    end
    MSD_tracks{1,i}=meanmsd;
    
end
% Prepare for averaging
MSD_tracks(:,find(all(cellfun(@isempty, MSD_tracks),1))) = []; % delete all empty cells in MSD_tracks!
maxSize = max(cellfun(@numel, MSD_tracks));    %# Get the maximum vector size
fcn = @(x) [x  nan(1,maxSize-numel(x))];  %# semicolon here
rmat = cellfun(fcn, MSD_tracks,'UniformOutput',false);  %# Pad each cell with NaNs
rmat=rmat';
plotMSD=cell2mat(rmat);


MSD_to_mean=[];
for i=1:length(plotMSD)
    if plotMSD(i,1)~= 0 && plotMSD(i,1)>0.001
        MSD_to_mean=[MSD_to_mean; plotMSD(i,:)];
    end
end

[Anzahltracks,bla]=size(MSD_to_mean);


f1=figure();
figure(f1);
linspace=[0:maxSize-1]*DeltaT;
for i=1:Anzahltracks
    loglog(linspace(1:end),MSD_to_mean(i,1:end),'color',[0.8 0.8 0.8]);
    hold on;
end

%Average
MSD_alldata=nanmean(MSD_to_mean(:,1:end));
MSD_weights=1./(nanstd(MSD_to_mean(:,1:end)));

figure(f1)
loglog(linspace(1:end),MSD_alldata(1:end),'r');
hold off


%% fit to random walk model to get persistence time and velocity!

% Adjust which span of your data you would like to analyse
[xData, yData, weights] = prepareCurveData(linspace(1:50),MSD_alldata(1:50), MSD_weights(1:50)); %here it is for the first 50 frames which corrpesonds to the first 5 s
% Set up fittype and options.
ft = fittype( '2*a*b^2*(x-a*(1-exp(-(x/a))))+c ', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Algorithm','Levenberg-Marquardt' );
opts.Display = 'Off';
opts.StartPoint = [1.3, 1.0 0.01]; % to improve fitting
opts.Weights = weights;

% Fit model to data5.
[fitresult, gof] = fit( xData, yData, ft, opts );

coefs=coeffvalues(fitresult);
% Plot fit with data.

f2=figure( 'Name', ' Fit' );
figure(f2);
h =plot( fitresult, xData, yData);
hold off

f3=figure('Name','LogFit');
figure(f3);
for i=1:Anzahltracks
    plot(linspace(1:200),MSD_to_mean(i,1:200),'color',[0.8 0.8 0.8]);
    hold on;
end
errorbar(linspace(1:200),MSD_alldata(1:200),nanstd(MSD_to_mean(:,1:200)),'-k');
hold on
xd=linspace(1:200);
plot(xd,(2*coefs(1)*coefs(2)^2*(xd-coefs(1)*(1-exp(-(xd/coefs(1)))))+coefs(3)),'color','#A2142F','Linewidth',2);
% plot(xd,(2*coefs(1)*coefs(2)^2*(xd-coefs(1)*(1-exp(-(xd/coefs(1)))))),'color','#A2142F','Linewidth',2);
hold on;
xd=linspace(1:200);
% plot(xd,(2*1.3*1.2^2*(xd-1.3*(1-exp(-(xd/1.3))))+0.01),'color','g','Linewidth',2);
set(gca, 'XScale','log', 'YScale','log')
set(gca,...
    'FontSize',10,...
    'FontWeight','bold',...
    'FontName','Arial')
ax=gca;
ax.LineWidth=1.5


tau=coefs(1)
v=coefs(2)
%% also for the single tracks, so you can choose

figure();
MSDplots=[];
%fit to MSD funktion to get correlation time for each single measurement
for i=1:size(MSD_to_mean,1)
    
    [xData, yData] = prepareCurveData(linspace(1:50),MSD_to_mean(i,1:50)); % adjust when needed
    
    % Set up fittype and options.
    ft = fittype( '2*a*b^2*(x-a*(1-exp(-(x/a))))+c ', 'independent', 'x', 'dependent', 'y' );
    
    opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Algorithm','Levenberg-Marquardt' );
    opts.Display = 'Off';
    opts.StartPoint = [1.3, 1.0 0.01];
    opts.Lower = [0.2, 0.8, 0];
    opts.Upper = [4, 3, 0.2];
    % opts.Weights = weights;
    
    % Fit model to data5.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    
    MSDtracksused{i}=yData;
    
    coefs=coeffvalues(fitresult);
    alpha = 0.95;
    ci = confint(fitresult);
    
    se = (ci(2,:)-ci(1,:))/2 % Standard Error
    
    
    if se(1)<coefs(1)*0.5 && se(2)<coefs(2)*0.5
        % Plot fit with data.
        
        
        %  h =plot( fitresult, xData, yData);
        % hold on
        
        MSDplots=[MSDplots;MSD_to_mean(i,1:50)]
        
        taui{i}=coefs(1);
        vi{i}=coefs(2);
        a{i}=coefs(3);
    end
end
hold off

ai=cell2mat(a);
meanai=mean(ai)
vi=cell2mat(vi)
meanv=mean(vi)
stdv=std(vi)/sqrt(length(vi))

figure();
hist(vi)

taui=cell2mat(taui)
meantau=mean(taui)
stdtau=std(taui)/sqrt(length(taui))

figure();
hist(taui)




f3=figure('Name','LogFit');
figure(f3);
for i=1:length(taui)
    loglog(linspace(1:50),MSDplots,'color',[0.8 0.8 0.8]);
    hold on;
end
hold on
xd=linspace(1:50);
loglog(xd,(2*meantau*meanv^2*(xd-meantau*(1-exp(-(xd/meantau))))+meanai),'color','#A2142F','Linewidth',2);
hold on;
set(gca, 'XScale','log', 'YScale','log')
set(gca,...
    'FontSize',10,...
    'FontWeight','bold',...
    'FontName','Arial')
ax=gca;
xlim([0.1 5]);
ax.LineWidth=1.5






% final plot
f9=figure();

figure(f9);
set(gcf,'position',[10,10,400,200])
xd=linspace(1:50);

loglog(xd,MSDplots,'color',[0.8 0.8 0.8]);
hold on;

xd=linspace(1:50);
loglog(xd,(2*meantau*meanv^2*(xd-meantau*(1-exp(-(xd/meantau))))+meanai),'color','#A2142F','Linewidth',2);
% plot(xd,(2*coefs(1)*coefs(2)^2*(xd-coefs(1)*(1-exp(-(xd/coefs(1)))))),'color','#A2142F','Linewidth',2);
hold on;

% plot(xd,(2*1.3*1.2^2*(xd-1.3*(1-exp(-(xd/1.3))))+0.01),'color','g','Linewidth',2);
set(gca, 'XScale','log', 'YScale','log')
set(gca,...
    'FontSize',10,...
    'FontWeight','bold',...
    'FontName','Arial')
xlabel(['Time [s]'])
ylabel(['MSD [µm^2]'])
legend(['Data'], ['Fit'],'Location', 'northwest')
ax=gca;
set(gca,'Box','on');
set(gca, 'layer', 'top');
xlim([0.099 5]);
ax.LineWidth=1.5

