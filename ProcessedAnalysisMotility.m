%% Average data of motility analysis
clear all;
clc;
close all;

cmap=[[128 128 128]./256; [70 70 70]./256; [254 128 0]./256;[254 78 0]./256];

%% Read Data

datafiles = dir(fullfile('yourFolder\*.mat'))



for ii = 1:length(datafiles)
    tmp = load(fullfile(datafiles(ii).folder, datafiles(ii).name));
    
    datastr = sprintf('Auswertung_%u', ii);  % Generate data string
    data.(datastr) = tmp
end

f1=figure();

xlabel('Radii from center [µm]');
ylabel('Diffusion constant [µm^2/s]');

B=cell2mat(struct2cell(data))
helpdata=data;
MeanMatrix=zeros(25,size(18,2));
% Bring data in right format
for j=1:length(datafiles)
    
    Xhelp=flip(B(j).XA);
    
    Y=flip(B(j).YA);
    Xmax=max(Xhelp);
    X=[];
    for i=1:length(Xhelp)
        X=[X Xmax-Xhelp(i)+1]
    end
    
    YA=B(j).ErrorYA;
    
    plot(X-1,Y,'color',[.5 .5 .5]);
    
    hold on;
    MeanMatrix(j,X)=Y;
    
end
hold on

% Average your data:
MeanMatrix(MeanMatrix==0)=NaN;
MeanalongMatrix=nanmean(MeanMatrix,1);
stdMatrix= nanstd(MeanMatrix,[],1);
temp=MeanalongMatrix;
temp(~isnan(temp)) = 1;
temp(isnan(temp)) = 0;
temp = find(temp);
first_non_NaN_index_of_X = temp(1);
linspace=[0:length(temp)-1]
Mean=MeanalongMatrix(~isnan(MeanalongMatrix));
STD=stdMatrix(~isnan(stdMatrix));
errorbar(temp-1, Mean,STD,'b');



%Plot the data
f3=figure('Renderer', 'painters', 'Position', [100 100 330 400]);
errorbar(temp-1, Mean,STD,'LineWidth',1.5,'Color',cmap(1,:),'Marker','x');
hold on
ylabel(['D [$\mu m^2 s^{-1}$]']);
xlabel(['d [$\mu m$]']);
xlim([-1 20])
set(gca,...
    'FontSize',10,...
    'FontWeight','bold',...
    'FontName','Arial')

ax=gca;
ax.LineWidth=1.5
legend('data')

