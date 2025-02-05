%% Program to evaluate Analysed Data from Interactionstates.m
% Example for four conditions:

clear all;
close all;
clc;
% number of colors in colormap depends on number of conditons:
cmap=[[160 160 160]./256; [27 120 55]./256;[90 174 97]./256;[172 211 158]./256];

set(0,'DefaultAxesColorOrder',cmap)
set(0,'defaultAxesFontSize',10);
set(groot, 'DefaultAxesFontName', 'Arial');

%% DoubleTrapAnalysiswithrightParameters

fname = {fullfile('Folderconditon1\','**', 'Auswertung.mat'),...
    fullfile('Folderconditon2\','**', 'Auswertung.mat'),...
    fullfile('Foldercondition3\','**', 'Auswertung.mat'),...
    fullfile('Folderconditon4\','**', 'Auswertung.mat')};

% Read in the data:
for co=1:length(fname)
    %% read files
    datafiles = dir(fname{co});
    
    for ii = 1:length(datafiles)
        tmp = load(fullfile(datafiles(ii).folder, datafiles(ii).name));
        
        datastr = sprintf('Auswertung_%u', ii);  % Generate data string
        data.(datastr) = tmp;
    end
    
    
    %% convert files
    B=cell2mat(struct2cell(data));
    %% Parameters you like to analyze
    
    Frequencieseventsall=[];
    %determine duration
    durationall=[];
    % determine frequencies of interaction states
    frequbundlall=[];
    frequpausall=[];
    frequretrall=[];
    frequeloall=[];
    %forces
    D=[];
    % how often appeared a state
    Leniret=[];
    Lenibundl=[];
    LeniPausi=[];
    LeniElod=[];
    Leniretd=[];
    Lenibundld=[];
    LeniPausid=[];
    LeniElo=[];
    % average force
    meanrupg4=[];
    %to plot the rupture force dustribution
    Forcehistogram=[];
    % all forces
    Dunsauber=[];
    % exclude forces larger than 80 pN
    Due=[];
    
    %% Forces
    
    for i=1:length(datafiles)
        C=B(i).Ruptureforces;
        C=cell2mat(C);
        if isempty(C)==1
            C=[1 100];
            
        end
        Dunsauber=[Dunsauber;C];
        for j=1:size(C,1)
            
            if C(j,2)<80 % Exclude forces beyond 80 pN
                Due=[Due; C(j,2)];
            end
            
            
            
            if C(j,2)>=80
                C(j,2)=100;
                D=[D; C(j,2)];
            elseif C(j,2)>1
                D=[D; C(j,2)];
                
            end
            
            
        end
        
        %% Probabilities and frequencies of interaction states
        
        
        Frequenciesevents=B(i).cuttedprobnoev;
        duration=B(i).durationcutted;
        % already normalized to whole duration of track
        frequbundl=B(i).cuttedprobbundl;
        frequpaus=B(i).cuttedprobpaus;
        frequretr=B(i).cuttedprobrup;
        frequelo=B(i).cuttedprobelo;
        
        
        if isempty(Frequenciesevents)==0
            if duration >30 && Frequenciesevents >0.1 && frequretr<0.3 % Exclude short tracks ov less than 30 s to analyse frquencies of events
                % to average over all tracks
                Frequencieseventsall=[Frequencieseventsall;1-Frequenciesevents]
                durationall=[durationall;duration];
                frequbundlall=[frequbundlall;frequbundl];
                frequpausall=[frequpausall;frequpaus];
                frequretrall=[frequretrall;frequretr];
                frequeloall=[ frequeloall;frequelo];
                if isempty(B(i).Rupturestate)==1
                    leni=1;
                else
                    
                    leni=length(B(i).Rupturestate);
                end
                lenieloi=length(B(i).Elongationstate);
                lenipause=length(B(i).Pausingstate);
                lenibundli=length(B(i).Bundlingstate);
                Leniret=[Leniret;leni]; % how often was retracted (count this array)
                
                Lenibundl=[Lenibundl;lenibundli];
                LeniPausi=[LeniPausi;lenipause];
                LeniElo=[LeniElo;lenieloi];
                
                
                Leniretd=[Leniretd;leni/duration]; % normalized to duration
                
                Lenibundld=[Lenibundld;lenibundli/duration];
                LeniPausid=[LeniPausid;lenipause/duration];
                LeniElod=[LeniElod;lenieloi/duration];
                
                
                
                
            end
        end
        
    end
    
    
    % average
    Freqret{co}=nanmean(Leniretd)
    Freqbundl{co}=nanmean(Lenibundld)
    FreqPausi{co}=nanmean(LeniPausid)
    FreqElo{co}=nanmean(LeniElod)
    
    
    erFreqret{co}=nanstd(Leniretd)/sqrt(length(~isnan(Leniretd)))
    erFreqbundl{co}=nanstd(Lenibundld)/sqrt(length(~isnan(Lenibundld)))
    erFreqPausi{co}=nanstd(LeniPausid)/sqrt(length(~isnan(LeniPausid)))
    erFreqElo{co}=nanstd(LeniElod)/sqrt(length(~isnan(LeniElod)))
    
    Probtobeboundall{co}=Frequencieseventsall;
    
    durationallover5{co}=durationall;
    
    durationsfornum{co}=sum(durationall);
    
    bundproball{co}=frequbundlall;
    pausingproball{co}=frequpausall;
    retraproball{co}=frequretrall;
    eloproball{co}=frequeloall;
    Ruptureforces{co}=D;
    Ruptureforceswithgoodinterval{co}=Due;
    Ruptureforcesraw{co}=Dunsauber;
    
    
    Nret{co}=sum(Leniret);
    
    Nbundl{co}=sum(Lenibundl);
    Npaus{co}=sum(LeniPausi);
    NElo{co}=sum(LeniElo);
    
    
    retrprob{co}=nanmean(frequretrall);
    eloprob{co}=nanmean(frequeloall);
    pausprob{co}=nanmean(frequpausall);
    bundlprob{co}=nanmean(frequbundlall);
    errretrprob{co}=nanstd(frequretrall)/sqrt(length(~isnan(frequretrall)));
    erreloprob{co}=nanstd(frequeloall)/sqrt(length(~isnan(frequeloall)));
    errpausprob{co}=nanstd(frequpausall)/sqrt(length(~isnan(frequpausall)));
    errbundlprob{co}=nanstd(frequbundlall)/sqrt(length(~isnan(frequbundlall)));
    errFrequenciesstates{co}=[errretrprob; erreloprob; errpausprob; errbundlprob];
    
    
    Frequenciesstates{co}=[retrprob; eloprob; pausprob; bundlprob];
    probtobebound{co}=nanmean(Frequencieseventsall);
    errproptobebound{co}=nanstd(Frequencieseventsall)/sqrt(length(~isnan(Frequencieseventsall)));
    
    
    
    
    Allforces{co}=D;
    fclose('all');
    clear B
    clear datafiles
    clear tmp
    clear datastr
    clear data
end
Legends={'a','b','c','d'}
%% Probability to be bound (excluded no event data)
f1=figure();
figure(f1);
edgesbound=[0:0.1:1];
for i=1:length(fname)
    
    if i==1
        histogram(Probtobeboundall{i},edgesbound,'Normalization','probability','EdgeColor','none','FaceColor',cmap(1,:),'FaceAlpha',1);
        hold on;
    elseif i==2
        histogram(Probtobeboundall{i},edgesbound,'Normalization','probability','DisplayStyle','stairs','LineStyle',':','Linewidth',2);
        hold on;
    else
        
    end
    
end
legend(Legends);
hold off
%% Probability that a specific state happened

Frequ=figure('Renderer', 'painters', 'Position', [100 100 330 400]);
figure(Frequ);
Freque=[];

A=Frequenciesstates{1,length(fname)};

A=cell2mat(A);
A=[A(:,1) A(:,2) A(:,3) A(:,4)];
bari=bar(A,0.9,'EdgeColor','none');
bari(1).FaceColor=cmap(1,:)
bari(2).FaceColor=cmap(2,:)
bari(3).FaceColor=cmap(3,:)
bari(4).FaceColor=cmap(4,:)

errFrequ= errFrequenciesstates{1,length(fname)};
errFrequ=cell2mat(errFrequ);
errFrequ=[errFrequ(:,1) errFrequ(:,2)  errFrequ(:,3) errFrequ(:,4)];
hold on;
ngroups = size(A, 1);
nbars = size(A, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, A(:,i), errFrequ(:,i), '.','Color',[0 0 0],'LineWidth',1.5);
end

% hold on;
% errorbar(Frequencies,errFrequencies)
% title('Frequencies of states')
ylabel(['p'])
xticklabels({'Retraction','Elongation','Pausing','Bundling'})
xtickangle(45)
legend(Legends);
set(gca,...
    'FontSize',10,...
    'FontWeight','bold',...
    'FontName','Arial')
% set(gca, 'Position', get(gca, 'OuterPosition') - ...
%     get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
ax=gca;
ax.LineWidth=1.5
hold off


hold off

%% Probability of frequencies of events per second

Frequs=figure('Renderer', 'painters', 'Position', [100 100 330 400]);;
figure(Frequs);
Freques=[];



As=[cell2mat(Freqret); cell2mat(FreqElo);cell2mat(FreqPausi);cell2mat(Freqbundl)]
As=[As(:,1) As(:,2) As(:,3) As(:,4)];
bari=bar(As,0.9,'EdgeColor','none');
bari(1).FaceColor=cmap(1,:)
bari(2).FaceColor=cmap(2,:)
bari(3).FaceColor=cmap(3,:)
bari(4).FaceColor=cmap(4,:)


errFrequs=[cell2mat(erFreqret); cell2mat(erFreqElo);cell2mat(erFreqPausi);cell2mat(erFreqbundl)]
errFrequs=[errFrequs(:,1) errFrequs(:,2)  errFrequs(:,3) errFrequs(:,4)];
hold on;
ngroups = size(As, 1);
nbars = size(As, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, As(:,i), errFrequs(:,i), '.','Color',[0 0 0],'LineWidth',1.5);
end

% hold on;
% errorbar(Frequencies,errFrequencies)
% title('Frequencies of states')
ylabel(['Events/s'])
xticklabels({'Retraction','Elongation','Pausing','Bundling'})
xtickangle(45)
legend(Legends);
set(gca,...
    'FontSize',10,...
    'FontWeight','bold',...
    'FontName','Arial')
% set(gca, 'Position', get(gca, 'OuterPosition') - ...
%     get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
ax=gca;
ax.LineWidth=1.5
hold off


hold off






%% Frequencies of events

Frequenciiiies=[];
Number=figure();
figure(Number);
Numbers=[];



for i=[1 2 3 4]
    
    Frequenciiiies{i}=[Nret{i}/durationsfornum{i};  Nbundl{i}/durationsfornum{i}; NElo{i}/durationsfornum{i}; Npaus{i}/durationsfornum{i}];
    
end


A=cell2mat(Frequenciiiies);

bari=bar(A);
bari(1).FaceColor=cmap(1,:)
bari(2).FaceColor=cmap(2,:)
bari(3).FaceColor=cmap(3,:)
bari(4).FaceColor=cmap(4,:)

xticklabels({'Retraction','Elongation','Pausing','Bundling'})
legend(Legends');
hold off





%% plot rupture force distribution (all data)
edges=[0:5:170];
f4=figure();
figure(f4);
for i=1:length(fname)
    if i==1
        histogram(Ruptureforcesraw{i},edges,'Normalization','probability','EdgeColor','none','FaceColor',cmap(1,:),'FaceAlpha',0.3);
        hold on;
    elseif i==2
        histogram(Ruptureforcesraw{i},edges,'Normalization','probability','DisplayStyle','stairs','LineStyle',':','Linewidth',2);
        hold on;
    else
        histogram(Ruptureforcesraw{i},edges,'Normalization','probability','DisplayStyle','stairs','LineStyle','-','Linewidth',2);
        hold on
        
        
    end
end
legend(Legends,'box', 'off');
hold off

edges=[0:5:110];

f5=figure('Renderer', 'painters', 'Position', [100 100 400 250]);


figure(f5);


for i=1:length(fname)
    if i==1
        histogram(Ruptureforces{i},edges,'Normalization','probability','EdgeColor','none','FaceColor',cmap(1,:),'FaceAlpha',1);
        hold on;
    elseif i==2
        histogram(Ruptureforces{i},edges,'Normalization','probability','Edgecolor','none','FaceColor',cmap(2,:),'FaceAlpha',1);
        hold on;
        
    else
        histogram(Ruptureforces{i},edges,'Normalization','probability','DisplayStyle','stairs','LineStyle','-','Linewidth',2);
        hold on;
        
        
    end
end
% title('Rupture forces')
xlabel(['Force [pN]'])
ylabel(['p'])
% legend(Legends,'box', 'on','LineWidth', 0.2,'FontWeight','normal');
legend('off')
ylim([0 0.9])
set(gca,...
    'FontSize',10,...
    'FontWeight','bold',...
    'FontName','Arial')
% set(gca, 'Position', get(gca, 'OuterPosition') - ...
%     get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
ax=gca;
ax.LineWidth=1.5
hold off
%% plot rupture force distribution (cut-off beyond 80 pN)
meanruptures=[];
errormeanruptureforces=[];
for i=1:length(fname)
    
    meanrupture=mean(Ruptureforceswithgoodinterval{i});
    stdmeanrupture=std(Ruptureforceswithgoodinterval{i})/sqrt(length(Ruptureforceswithgoodinterval{i}));
    meanruptures=[meanruptures; meanrupture];
    errormeanruptureforces=[errormeanruptureforces;stdmeanrupture];
end


figure('Position', [100 100 300 300]);

barpi=bar(meanruptures',0.9,'EdgeColor','k','FaceColor','flat');
barpi.CData(1,:)=cmap(1,:)
barpi.CData(2,:)=cmap(2,:)
barpi.CData(3,:)=cmap(3,:)
barpi.CData(4,:)=cmap(4,:)


hold on;
ngroups = size(meanruptures, 1);
nbars = size(meanruptures, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, meanruptures(:,i), errormeanruptureforces(:,i), '.','color', 'k');
end

ylabel('Rupture force [pN]')
xticklabels(Legends);
xtickangle(45)
set(gca,...
    'FontSize',10,...
    'FontWeight','bold',...
    'FontName','Arial')
% set(gca, 'Position', get(gca, 'OuterPosition') - ...
%     get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
ax=gca;
ax.LineWidth=1.5
hold off


%% mean forces whole interval


meanrupturesraw=[];
errormeanruptureforcesraw=[];
for i=1:length(fname)
    
    meanruptureraw=mean(Ruptureforcesraw{i});
    stdmeanruptureraw=std(Ruptureforcesraw{i});
    meanrupturesraw=[meanrupturesraw; meanruptureraw];
    errormeanruptureforcesraw=[errormeanruptureforcesraw;stdmeanruptureraw];
end

forcesraw=figure();
figure(forcesraw);
bar(meanrupturesraw(:,2));
meanrupturesraw=meanrupturesraw(:,2);
errormeanruptureforcesraw=errormeanruptureforcesraw(:,2);
hold on;
ngroups = size(meanrupturesraw, 1);
nbars = size(meanrupturesraw, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, meanrupturesraw(:,i), errormeanruptureforcesraw(:,i), '.','color', cmap(i,:));
end

title(['Forcesraw'])
xticklabels(Legends);
hold off



%% kstest für rupture forces:

[h1,p1]=kstest2(Ruptureforceswithgoodinterval{2},Ruptureforceswithgoodinterval{4})
