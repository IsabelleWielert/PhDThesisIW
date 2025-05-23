%%Script: Moving mean of diffusionconstants in colonies (radially resolved)
% You can use this script to just plot data or even to get the mean of data
% Name from developmental stage: meanoftimepointsfromshell

set(groot, 'DefaultTextInterpreter', 'tex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'tex');
set(groot, 'DefaultLegendInterpreter', 'tex');
% Schriftart und Größe global setzen
set(groot, 'DefaultAxesFontName', 'Arial');
set(groot, 'DefaultTextFontName', 'Arial');
set(groot, 'DefaultLegendFontName', 'Arial');
set(groot, 'DefaultAxesFontSize', 10);
set(groot, 'DefaultTextFontSize', 10);
set(groot, 'DefaultLegendFontSize', 10);
clear all;
clc;
close all
% general settings
set(groot, 'DefaultAxesFontName', 'Arial');
set(0,'DefaultAxesColorOrder',hsv(16));
cmap=[[216 27 96]./255;[30 136 229]./255; [255 193 7]./255; [0 77 64]./255; [254 97 0]./255 ;[120 94 240]./255;[120 94 100]./255;[120 94 10]./255];



% Hauptpfad setzen (anpassen!)
hauptordner = 'Yourfolder';

% Infos über alle Dateien und Unterordner holen
allePfadInfos = dir(fullfile(hauptordner, '**', '*','*')); % '**' = rekursiv

% Nur Verzeichnisse filtern (ohne '.' und '..')
alleOrdner = allePfadInfos([allePfadInfos.isdir]);
alleOrdner = alleOrdner(~ismember({alleOrdner.name}, {'.', '..'}));

Radius_innercore=[];
Widthshell=[];
foldinginnercore=[];
% Loop über alle Ordner

r30min=[];
r60min=[]; 
r90min=[];
for i = 1:length(alleOrdner)
    aktuellerOrdner = fullfile(alleOrdner(i).folder, alleOrdner(i).name);
    fprintf('Ordner: %s\n', aktuellerOrdner);
    
    % Beispiel: Alle .mat-Dateien in diesem Unterordner finden
    matFiles = dir(fullfile(aktuellerOrdner, '*.mat'));
    
      % Die .mat-Dateien nach Namen sortieren
          matFileNames = string({matFiles.name});
    [~, idx] = sort(matFileNames,'descend');  % Sortiere nach Namen
    matFiles = matFiles(idx);          % Wende die Sortierung an
    
    names={'30 min', '60 min', '90 min', '120 min'};
    for j = 1:length(matFiles)
        
 
        dateipfad = load(fullfile(matFiles(j).folder, matFiles(j).name));
        
    
        datastr = sprintf('Auswertung_%u', j);  % Generate data string
        data.(datastr) =dateipfad;
%        names{j}= matFiles(j).name
       
    end
    B=cell2mat(struct2cell(data));
    
    
    % Plot the data
    figure('Units', 'Inches', 'Position', [4,4, 4, 2.5]);
%0.1*length(B(k).newYA)
 for k=1:length(B)
    
%     plot(B(k).XA_sortedfromedge,max(smoothdata(B(k).newYA,'movmean',0.2*length(B(k).newYA)))-smoothdata(B(k).newYA,'movmean',0.2*length(B(k).newYA)), 'LineWidth',1.5,'Color', cmap(k,:));
%     hold on;
plot(max(B(k).XA_sortedfromedge)-B(k).XA_sortedfromedge,smoothdata(B(k).newYA,'movmean',0.2*length(B(k).newYA)), 'LineWidth',1.5,'Color', cmap(k,:));
hold on; 
%     legend('off');
    ylabel('D [µm^2/s]');
    xlabel('Distance from centre of colony [µm]');
ylim([0, 0.0015])
xlim([0, 25])

%      plot(max(B(k).XA_sortedfromedge)-B(k).XA_sortedfromedge,smoothdata(B(k).newYA,'movmean',0.2*length(B(k).newYA)), 'LineWidth',1.5,'Color', cmap(k,:));
      legend(names,'Box', 'off', 'Location', 'northeast')
legend off
end
presortforplot=[];
for k=1:length(B)
       % radius of the motile core
     figure('Units', 'Inches', 'Position', [7, 7, 4, 2.5]);
    findpeaks(max(smoothdata(B(k).newYA,'movmean',0.2*length(B(k).newYA)))-smoothdata(B(k).newYA,'movmean',0.2*length(B(k).newYA)),B(k).XA_sortedfromedge,'MinPeakProminence',0.0002,'Annotate','extents')
    hold on
    grid off
        ylabel('D_{max}-D [µm^2/s]');
    xlabel('Distance from edge of colony [µm]');
    hold off
    [pks,locs,widths,proms]=findpeaks(max(smoothdata(B(k).newYA,'movmean',0.2*length(B(k).newYA)))-smoothdata(B(k).newYA,'movmean',0.2*length(B(k).newYA)),B(k).XA_sortedfromedge,'MinPeakProminence',0.0002,'Annotate','extents');
    if ~isempty(pks)
    ax = gca;
    lines = ax.Children;
    x = lines(1).XData';
    x = x(~isnan(x));
    x = transpose(reshape(x,2,[]));
    
    
    % Determine Radius of inner motile core
    Radius_innercore=[Radius_innercore; max(B(k).XA_sortedfromedge)-max(x, [], 'all')]
%     Radius_innercore=[Radius_innercore;min(x, [], 'all')]
    Widthshell=[Widthshell; max(widths)]
%  Widthshell=[Widthshell; max(x, [], 'all')]
%     presortforplot=[presortforplot; max(B(k).XA_sortedfromedge)-max(x, [], 'all'),max(widths)]
 presortforplot=[presortforplot; max(B(k).XA_sortedfromedge)-max(x, [], 'all'),max(widths)]
 
if k==1
    r30min=[r30min; max(B(k).XA_sortedfromedge)-max(x, [], 'all') ]
elseif k==2
   r60min=[r60min; max(B(k).XA_sortedfromedge)-max(x, [], 'all') ]
elseif k==3
       r90min=[r90min; max(B(k).XA_sortedfromedge)-max(x, [], 'all') ]
end
    end
  
end
[val,indexi]=max(presortforplot)
foldinginnercore=[foldinginnercore; presortforplot(indexi(1),:)]


hold off;
clear data
    
end
%%


figure(); 

scatter(foldinginnercore(1:end-1,1),foldinginnercore(1:end-1,1)+foldinginnercore(1:end-1,2))

%% read files
% datafiles = dir(fname{co});
% 
% for ii = 1:length(datafiles)
%     tmp = load(fullfile(datafiles(ii).folder, datafiles(ii).name));
%     
%     datastr = sprintf('Auswertung_%u', ii);  % Generate data string
%     data.(datastr) = tmp;
% end
% 
% 
% % read out the data
% B=cell2mat(struct2cell(data));
% 
% 
% 
% 
% % Plot the data
% figure('Position', [300 400 600 400]);
% 
% for i=1:length(datafiles)
%     
%     plot(B(i).XAsortedfromedge,smoothdata(B(i).newYA,'movmean',0.1*length(B(i).newYA)), 'LineWidth','color', cmap(i,:));
%     hold on;
%     legend('Location', 'northwest');
%     ylabel('Diffusion constant [µm^2/s]');
%     xlabel('Distance from edge of colony [µm]');
%     title('Radial motility (X hours prior to folding)');
% 
%     
% end
% hold off;
% 
% 
% 
% Radius_innercore=[];
% Widthshell=[];
% 
% for i=1:length(datafiles)
%     
%     % findpeaks and check in plot if this is working, otherwise adjust Peak
%     % Prominence. The data is inversed and you can determine a peak of low
%     % motility.Then, substracting half-width from radius, resulting in inner
%     % radius of the motile core
%     figure();
%     findpeaks(max(smoothdata(B(i).newYA,'movmean',0.1*length(B(i).newYA)))-smoothdata(B(i).newYA,'movmean',0.1*length(B(i).newYA)),B(i).XA_sorted,'MinPeakProminence',0.0002,'Annotate','extents')
%     hold off
%     [pks,locs,widths,proms]=findpeaks(max(smoothdata(B(i).newYA,'movmean',0.1*length(B(i).newYA)))-smoothdata(B(i).newYA,'movmean',0.1*length(B(i).newYA)),B(i).XA_sorted,'MinPeakProminence',0.0002,'Annotate','extents');
%     
%     ax = gca;
%     lines = ax.Children;
%     x = lines(1).XData';
%     x = x(~isnan(x));
%     x = transpose(reshape(x,2,[]));
%     
%     
%     % Determine Radius of inner motile core
%     %Radius_innercore=[Radius_innercore; max(B(i).XA_sortedfromedge)-max(widths)]
%     Radius_innercore=[Radius_innercore;min(x, [], 'all')]
%     Widthshell=[Widthshell; max(widths)]
%     
%     
% end
% 
% 
  figure('Units', 'Inches', 'Position', [4,4, 4, 2.5]);
daten1=r90min'
daten2=r60min'
daten3=r30min'

alleDaten = [daten1, daten2, daten3];
gruppen = [ones(1,length(daten1)), ...
           2*ones(1,length(daten2)), ...
           3*ones(1,length(daten3))];

% Boxplot zeichnen
boxplot(alleDaten, gruppen, 'Labels', {'90 min', '60 min', '30 min'}, 'MedianStyle','line', 'Colors','k')
hold on



% Farben für Punkte
farben = [[216 27 96]./255;[30 136 229]./255; [255 193 7]./255];

%
% Rohdaten als farbige Punkte plotten
rng(1); % für Reproduzierbarkeit
x1 = 1 + 0.1*randn(1,length(daten1));
x2 = 2 + 0.1*randn(1,length(daten2));
x3 = 3 + 0.1*randn(1,length(daten3));

scatter(x1, daten1, 40, farben(1,:), 'filled', 'MarkerFaceAlpha', 0.7)
scatter(x2, daten2, 40, farben(2,:), 'filled', 'MarkerFaceAlpha', 0.7)
scatter(x3, daten3, 40, farben(3,:), 'filled', 'MarkerFaceAlpha', 0.7)
ylabel('<r> [µm]')


hold off
