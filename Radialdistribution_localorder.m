%% Program to calculate radial distribution function from single-cell
% tracking 
clear all;
clc;
cmap={[128 128 128]./256; [179 47 109]./256; [254 128 0]./256;[20 149 0]./256};

% Enter your folder name, where you saved the positions
allpath = {'Your folder names'};
% enter legend names
leg = {'x','x'};

% Depends on which conditions you investigated, needs to be the same size
% as your number of conditions
normdiameter = [0.51*2];

figure();
for u = 1:length(allpath)
path = allpath{u};

name = split(path,'\');
% Get a list of all files and folders in this folder.
files = dir(path);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories
subFolders = files(dirFlags);
subFolders(1:2) = [];

rsl = [2,5];
xyscale = 0.0792354/rsl(1);
zscale = 0.2/rsl(2);
distcell = {};
binsize = 0.075;
maxbin = 30;
edges = [0:binsize:maxbin];
radius = [0+binsize/2:binsize:maxbin-binsize/2];
% Nall = {};
% RDFcell = {};
ind = 1;
zthresh = 2;
h = 10-zthresh;
ind2 = 1;
Nnorm = 1000;
% dcell = {};
Nall = nan(length(radius),length(subFolders));
% figure();
% hold all;
for i = 1:length(subFolders)
    if exist([path,'\',subFolders(i).name,'\data.mat'], 'file') == 2

        load([path,'\',subFolders(i).name,'\data.mat'],'r');
        rgreen =  r{1,1};

        x = double(rgreen(:,1)).*xyscale;
        y = double(rgreen(:,2)).*xyscale;
        z = double(rgreen(:,3)).*zscale;
        [zsort,I] = sort(z);
        xsort = x(I);
        ysort = y(I);
        xcent = sum(x)/length(x);
        ycent = sum(y)/length(y);
        distc = zeros(1,length(x));
        for j = 1:length(x)
           distc(j) = sqrt((x(j)-xcent)^2+(y(j)-ycent)^2);
        end
        
        [n,edges2] = histcounts(distc,'BinWidth',0.25);
        [~,I] = max(n);
        rmax = edges2(I);
        if rmax > 2
        for j = 1:length(x)-1
           for k = j+1:length(x)
               dist1 = sqrt((x(j)-xcent)^2+(y(j)-ycent)^2);
               dist2 = sqrt((x(k)-xcent)^2+(y(k)-ycent)^2);
               if z(k) > zthresh && z(j) > zthresh 
                   tempdist = sqrt((x(k)-x(j))^2+(y(k)-y(j))^2+(z(k)-z(j))^2);
%                    if tempdist<rmax
                       dist(ind) = tempdist;
                       ind = ind+1;
%                    end
               end
           end
        end

clear x y z;
        
        N = histcounts(dist,edges);
        n_norm = movmean(N,15);


        Nall(:,i) = N./n_norm; 
%      
        end
        dist = [];
        ind = 1;
        clear s;
        clear x y z;
    end
end


RDF = nanmean(Nall,2);
RDFerr = nanstd(Nall,0,2);
numerr = sum(~isnan(Nall),2);

hold on;
e1=errorbar(radius(1:50),RDF(1:50),RDFerr(1:50)./sqrt(numerr(1:50)),'LineWidth',1.5);
    set(e1,'Color',cmap{u})
    set(e1,'MarkerEdgeColor',cmap{u})



end

xlabel('r [µm]');

ylabel('g(r)');
legend(leg);

        hold off;
        




% save your data
% save([path,'\RDF_Data.mat'],'distcell');