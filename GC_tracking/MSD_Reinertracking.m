%%%% Skript um MSD und correltaion time zu berechnen 
scale=0.83;%micron
DeltaT=1/10; %weil mit 8 Hz gemessen wurde 
MSD_to_mean=[];
files=dir('C:\Users\Isabelle\Documents\Twitching Assay\Data_all\G4\*.txt');
names={files.name};
for i=1:length(names)
%data=readtable('C:\Users\Isabelle\Documents\Twitching Assay\Daten\19_01_20\track_NG36.csv', 'HeaderLines',1);
fileID = fopen(['C:\Users\Isabelle\Documents\Twitching Assay\Data_all\G4\', names{i}],'r');
formatSpec = '%f %f %f %f %f %f';
sizeA = [5 Inf];
A = fscanf(fileID,formatSpec,sizeA);
datasorth=A';
%f
%sidn positionen schon in um ?


% erstmal nach track sortieren 

% 
% 
% datasort=sortrows(data,3);
% A=datasort{:,3};
% ID=str2double(A);
% 
% 
% datasorth=table2array(datasort(:,5:6));

PosID=[datasorth(:,3) datasorth(:,4) datasorth(:,1)];


%%
%Calculation of MSD 

  for i=1:max(PosID(:,3))-1

    px=PosID(PosID(:,3)==i,1);
    py=PosID(PosID(:,3)==i,2);
    meanmsd=zeros(1,length(px)+1);
    N=length(px);
    
    
    if N > 20
         %if N<150

        for j=1:N


            for k=j+1:N
                    dist=abs((px(j)-px(k))^2+(py(j)-py(k))^2);
                    h = k-j;
                    meanmsd(h)=meanmsd(h) + dist;

            end
             meanmsd(j)=meanmsd(j)/(N-j+1);

        end
             MSD_tracks{1,i}=meanmsd;
 
    %end

       
    end
  end
 
 MSD_tracks(:,find(all(cellfun(@isempty, MSD_tracks),1))) = []; % delete all empty cells in MSD_tracks!
 %MSD_tracks=cell2table(MSD_tracks);
 
%  MSD_tracks=vertcat(MSD_tracks{1,:});
maxSize = max(cellfun(@numel, MSD_tracks));    %# Get the maximum vector size
fcn = @(x) [x  nan(1,maxSize-numel(x))];  %# semicolon here
rmat = cellfun(fcn, MSD_tracks,'UniformOutput',false);  %# Pad each cell with NaNs
rmat=rmat';
plotMSD=cell2mat(rmat);

plotMSD=nanmean(plotMSD(:,1:end));
for i=1:1300
plotMSD=[plotMSD nan];
if length(plotMSD)==1200
    break
end
end
MSD_to_mean=[MSD_to_mean; plotMSD];
end

f1=figure(); 
figure(f1);
linspace=[1:1200]*DeltaT;
for i=1:length(names)
    
loglog(linspace(1:end),MSD_to_mean(i,1:end),'*r');
hold on;
end

MSD_alldata=nanmean(MSD_to_mean(:,1:end));
figure(f1)
loglog(linspace(1:end),MSD_alldata(1:end),'b');
hold off
% 
   
%fit to MSD funktion to get correlation time !


[xData, yData] = prepareCurveData(linspace(1:round(8*10)),MSD_alldata(1:8*10));

% Set up fittype and options.
ft = fittype( '2*a*b^2*(x-a*(1-exp(-(x/a)))) + c', 'independent', 'x', 'dependent', 'y' );

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [1,10,90];
opts.Lower = [0.1,2,10];
opts.Upper = [2,30,300];
% Fit model to data5.
[fitresult, gof] = fit( xData, yData, ft, opts );
coefs=coeffvalues(fitresult);
% Plot fit with data.

f2=figure( 'Name', ' Fit' );
figure(f2);
h =plot( fitresult, xData, yData);
hold off

f7=figure('Name','LogFit');
figure(f7);
loglog(xData,yData,'.')
hold on 
xd=linspace(1:80);
loglog(xd,(2*coefs(1)*coefs(2)^2*(xd-coefs(1)*(1-exp(-(xd/coefs(1)))))));
hold off

coefs(1)
coefs(2)
coefs(3)