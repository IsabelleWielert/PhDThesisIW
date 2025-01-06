%%Script: Moving mean of diffusionconstants in colonies (radially resolved)
% You can use this script to just plot data or even to get the mean of data
% Name from developmental stage: meanoftimepointsfromshell 


clear all;
clc;
close all;
% general settings
set(groot, 'DefaultAxesFontName', 'Arial');
set(0,'DefaultAxesColorOrder',hsv(16));
cmap=[[216 27 96]./255;[30 136 229]./255; [255 193 7]./255; [0 77 64]./255; [254 97 0]./255 ;[120 94 240]./255;[120 94 100]./255;[120 94 10]./255];    

% open file 
fname = {fullfile('FILE\*.mat')};
arr=[7,5,3]; % if you woild like to normalize the diffusion constant, this will get you the minimum



 for co=1:length(fname)
   %% read files
    datafiles = dir(fname{co});
    
    for ii = 1:length(datafiles)
        tmp = load(fullfile(datafiles(ii).folder, datafiles(ii).name));
        
        datastr = sprintf('Auswertung_%u', ii);  % Generate data string
        data.(datastr) = tmp;
    end
    
    
% read out the data 
    B=cell2mat(struct2cell(data));
    % preparation of arrays you would like to average
    ALLX=[];
    ALLY=[];
    ALLERR=[];
    
    % Plot the data 
    figure('Position', [300 400 600 400]);
  
for i=1:length(datafiles)
    
    % the next lines determines the normalization for you data, around an
    % value, which you set in arr 
     closest = interp1(B(i).XA_sortedfromedge,B(i).XA_sortedfromedge,arr(co),'nearest');
     indi=find(closest==B(i).XA_sortedfromedge);
     HelpY=B(i).newYA; 
     norm=HelpY(indi);
     
     % preparation of your data to average (average means (all in one array and then use the movmean function))
     ALLX=[ALLX,B(i).XA_sortedfromedge];
     ALLY=[ALLY,B(i).newYA];%./norm   % you can decide wether you want to normalize the data or not
     ALLERR=[ALLERR,B(i).ErrorYA];%./norm


 
    plot(B(i).XA_sortedfromedge,smoothdata(B(i).newYA,'movmean',0.1*length(B(i).newYA)), 'LineWidth', 1.5); 
    hold on;
    legend('Location', 'northwest');
    ylabel('Diffusion constant [µm^2/s]');
    xlabel('Distance from edge of colony [µm]');
    title('Radial motility (X hours prior to folding)');

end
 hold off; 
 
 
 % bring the data in the right form, sorted 
[X_sorted, orderX] = sort(ALLX);
X_sortedfromedge=max(X_sorted)-X_sorted;
sortedY =ALLY(orderX);
sortederr=ALLERR(orderX);



figure(); 
 errorbar(X_sorted,smoothdata(sortedY,'movmean',0.05*length(sortedY)), sortederr,'color', cmap(co,:), 'LineWidth', 1.5); 
%  hold on;
%  yline(min(smoothdata(sortedY,'movmean',0.1*length(sortedY))))
 hold on 
   legend('Location', 'northwest');
   ylabel('Mean Diffusion constant [µm^2/s]');
   xlabel('Distance from edge of colony [µm]');
   title('Radial motility (1.5 h prior to folding)');
   hold off;
   
    figure();
 plot(X_sorted,smoothdata(sortedY,'movmean',0.1*length(sortedY)), 'color', cmap(co,:), 'LineWidth', 1.5); 
  hold on;
%  yline(min(smoothdata(sortedY,'movmean',0.1*length(sortedY))))
%  hold on;
 
   legend({'Folding', 'minimum','Non-folding', 'minimum'},'Location', 'northeast');
   ylabel('Mean Diffusion constant [µm^2/s]');
   xlabel('Distance from edge of colony [µm]');
%    title('Radial motility (2 hours prior to folding)');
   hold off;
end
