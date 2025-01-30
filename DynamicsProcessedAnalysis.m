%% Programm to get elongation velocities & retraction velocities from distribution of retraction aand elongation velocities
% fit a normal distribution to the data and get the mean
clc;
clear;
close all;

%% Read out the data
Lengthdistri=[];
Dauerdistri=[];
Elongationdistri=[];
Retractiondistri=[];

datafiles{1,1}= dir(fullfile('your folder with **\.mat'));
for z=1:size(datafiles,2)
    datafilesi=datafiles{1,z} % If you have more than condition to compare, then you need to add also these pathways (add a new datafiles{x,x})
    for ii = 1:length(datafilesi)
        
        tmp = load(fullfile(datafilesi(ii).folder, datafilesi(ii).name));
        
        datastr = sprintf('Auswertung_%u', ii);  % Generate data string
        data.(datastr) = tmp;
    end
    
    B=cell2mat(struct2cell(data));
    %% Read out the duration of the pilus track

  Lebensdauer=[];
    for i= 1:length(B)
        AnzahlP=cell2mat(B(i).Lengths)
        if length(AnzahlP)==1
            Dauer=B(i).durationcutted;
            Lebensdauer=[Lebensdauer;Dauer];
        end
    end
    %% Read out and collect length from maxima of pilus tip tracks
        Piluslengths=[];
    for i=1:length(B)
        A=B(i).Lengths;
        Piluslengths=[Piluslengths;A];
    end
    LengthPili=cell2mat(Piluslengths);
    Lengthss=[];
    for i=1:length(LengthPili)
        if LengthPili(i)>0.7 %only count when larger than 0.7 due to imaging quality
            Lengthss=[Lengthss;LengthPili(i)];
        end
    end
    % get histogram (to fit exponential distribution, if you want)
    [h,cent]=hist(Lengthss,12);
    [N,edges]=histcounts(LengthPili,12);
    Bins=[];
    for i=1:length(edges)-1
        binni=edges(i+1)+edges(i)/2;
        Bins=[Bins;binni];
    end
    
    %% Elongation velocities
    
    Elongationdata=[];
    for i=1:length(B)
        A=B(i).Elongationstate;
        Elongationdata=[Elongationdata;A];
    end
    Elongation=cell2mat(Elongationdata);
    
    Elongation = Elongation(~isnan(Elongation(:,3)),3);
    
    
    %% Retraction velocities
    
    B=cell2mat(struct2cell(data));
    Retractiondata=[];
    for i=1:length(B)
        C=B(i).Retractionstate;
        Retractiondata=[Retractiondata;C];
    end
    retraction=cell2mat(Retractiondata);
    retraction = retraction(~isnan(retraction(:,3)),3);
    
    
    
    %% Fits via fitdist
    
    
    %Lengths
    Lengthdistri{z}=LengthPili;
    %Elongation velocities:
    meanelo=fitdist(Elongation,'Normal');
    Elongationdistri{z}=Elongation;
    Eloparam(z,1:2)=meanelo.Params(1,1:2);
    Confiint_elo=paramci(meanelo);
    Confielo(z,1:2)=Confiint_elo(1,1:2);
    
    %Retraction velocities:
    meanretr=fitdist(retraction,'Normal');
    Retractiondistri{z}=retraction;
    Retrparam(z,1:2)=meanretr.Params;
    Confiint_retr=paramci(meanretr);
    Confiretr(z,1:2)=Confiint_retr(1,1:2);
    
    % Duration of tracks:
    % meandauer=fitdist(Lebensdauer,'Normal');
    Dauerdistri{z}=Lebensdauer;
    % Dauerparam(z,1:2)=meandauer.Params;
    % Confiint_dauer=paramci(meandauer);
    % Confidauer(z,1:2)=Confiint_dauer(1,1:2);
    
    
    
end

% Plots to visualize the data
figure();hist(cell2mat(Elongationdata(:,3))); hold off
figure();hist(cell2mat(Retractiondata(:,3))); hold off

