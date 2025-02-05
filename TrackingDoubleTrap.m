
%% Program to track spherical particle pairs which are horizontially or vertically aligned (used to track trapped gonococci)

close all;
clc;


% uncomment if you want to check calibration with equipartition theorem (read out the stiffness in KX1, KX2,KY1,KY2)
% KX1=[];
% KY1=[];
% KX2=[];
% KY2=[];

% open your file: 
allpath='Path of folder where your experiment (.tif format)';
liste=dir(allpath);

% read out all folders with images
Paths={liste.name};
Paths=Paths(:,3:end);
B={};
for o=1:length(Paths)

pfadnamen=fullfile(allpath,Paths{o})
B{o}=pfadnamen;
end




for o=1:length(B)

path=B{o};
 
folderraw = path;
folderrawAuswertung = [path,'_Auswertung'];


slpos = strfind(folderraw,'\'); 
outname = [folderraw(slpos(end):end),'_rechts.txt'];
FileNames = dir(folderraw);
slpos2 = strfind(folderraw,'\'); 
outname2 = [folderraw(slpos2(end):end),'_links.txt'];
FileNames2 = dir(folderraw);
if exist(folderrawAuswertung, 'dir') == 0
      mkdir(folderrawAuswertung)
      
      
% remove dirs from Filenames
for i = length(FileNames):-1:1
    if FileNames(i).isdir == 1
        FileNames(i) = [];
    else
        [~,~,ext] = fileparts(FileNames(i).name);
        if ~strcmpi(ext, '.tiff')
            FileNames(i) = [];
        end
    end
end

%get Filename cell from structure for sorting
Files = cell(length(FileNames),1);
for i = 1:length(FileNames)
    Files{i} = FileNames(i).name;
end  
Files = sort_nat(Files);

% read Movie
temp = imread(fullfile(folderraw,Files{1}));

Movieee = zeros(size(temp,1),size(temp,2),length(FileNames),'uint8');
for i = 1:length(FileNames)
    temp = imgaussfilt(imread(fullfile(folderraw,Files{i})),2);
   Movieee(:,:,i) = squeeze(temp(:,:,1));
end



% Get initial bead position (you can check once if you recalibrated the trap, where the equilibrium positions are roughly)
%  f1=figure('Name','Where is the Bead?'); imagesc(sum(Movieee(:,:,:),3)); axis equal; colormap('gray')
% coord=round([100 100]);
coord = round([100 110]); %(probably) more precise window

% define roi where bead positions are calculated
border=3;
xs = max([coord(1)-40 1]);
xe = min([coord(1)+40-border size(Movieee(:,:,1),2)-border]);
ys = max([coord(2)-40 1]);
ye = min([coord(2)+40-border size(Movieee(:,:,2),1)-border]);

% define filter
fltr4img = [1 1 1 1 1; 1 2 2 2 1; 1 2 4 2 1; 1 2 2 2 1; 1 1 1 1 1];
fltr4img = fltr4img / sum(fltr4img(:));

 
% set parameters for hough grid
min_=2;
max_=12; % vorher 12 
grdthres=0.1;
fltr4LM_R=9;
multirad=0.8;


% perform tracking
%  figure(); % uncomment for live tracking
coords=zeros(size(Movieee,3),2);
coords2=zeros(size(Movieee,3),2);
for i = 1:size(Movieee,3)
    rawimg = Movieee(ys-border:ye+border,xs-border:xe+border,i);

    rawimg=imgaussfilt(rawimg);
    imgfltrd = filter2( fltr4img , rawimg );
    imgfltrd = imgfltrd(1+border:end-border,1+border:end-border);
    [~, circen, ~] = CircularHough_Grd(imgfltrd, [min_ max_], grdthres, fltr4LM_R, multirad); % applies Hough Grid tracking (see function)
    if isempty(circen)
        circen=nan(1,2);
      fprintf('x');
        

     elseif numel(circen)<3 
        circen=nan(1,2);
        fprintf('x');
    else
        
        %% uncomment if you have another alignment of the trapped particles
%         if circen(1,2)>circen(2,2)    % vertical alignment
%         coords(i,:)=circen(1,:);
%         coords2(i,:)=circen(2,:);
%         else circen(1,2)<circen(2,2)
%     
%           coords(i,:)=circen(2,:);
%         coords2(i,:)=circen(1,:);
        if circen(1,1)>circen(2,1)      % horizontal alignment
        coords(i,:)=circen(1,:);
        coords2(i,:)=circen(2,:);
        elseif circen(1,1)<circen(2,1)
    
          coords(i,:)=circen(2,:);
        coords2(i,:)=circen(1,:);

        end
    end
%%    uncomment if you would like to see the live tracking, but it slows down the tracking
%              imagesc(imgfltrd); 
%      hold on
%       plot(coords(i,1),coords(i,2),'+r');
%       hold on; 
%        plot(coords2(i,1),coords2(i,2),'+r');
%      hold off
%       drawnow
    
%%    
    
end


% remove wrong positions
ind = find(coords(:,1)==0);
coords(ind,:) = [];


ind2 = find(coords2(:,1)==0);
coords2(ind2,:) = [];




end


fid=fopen([folderrawAuswertung, outname],'at');
for k= 1:size(Movieee,3)-length(ind)
    
    fprintf(fid,'%s \n',[num2str(coords(k,:))]); 
    
end
fclose(fid);


fid2=fopen([folderrawAuswertung, outname2],'at');
for k= 1:size(Movieee,3)-length(ind2)
    
    fprintf(fid2,'%s \n',[num2str(coords2(k,:))]); 
    
end
fclose(fid2);


% test figure to check the track

figure('Position',[300 ,300, 500, 500],'name',B{o});
plot(coords(:,2)-mean(coords(:,2)),'.');
grid on;
axis square;
hold on;

%% figure if you would like to see the positon distribution of the top bead via a heatmap:
% figure('Position',[300 ,300, 500, 500],'name','HoughgridobererBeadheatmap');
% x=coords(:,2)-mean(coords(:,2));
% y=coords(:,1)-mean(coords(:,1));
% pts = linspace(-1.5,1.5, 30);
% N = histcounts2(y(:), x(:), pts, pts);
% imagesc(pts, pts, N);
% axis equal;
% 
% hold off
% figure('Position',[300 ,300, 500, 500],'name','HoughgriduntererBead');

plot(coords2(:,2)-mean(coords2(:,2)),'.');

hold off; 

%% figure if you would like to see the positon distribution of the bottom bead via a heatmap:

% figure('Position',[300 ,300, 500, 500],'name','HoughgriduntererBeadheatmap');
% x2=coords2(:,2)-mean(coords2(:,2));
% y2=coords2(:,1)-mean(coords2(:,1));
% pts = linspace(-1.5,1.5, 30);
% N2 = histcounts2(y2(:), x2(:), pts, pts);
% imagesc(pts, pts, N2); 
% axis equal;
% 

%% further evaluation (if you would like to test the equipartition theorem)

% resx = 88.652; 
% resy = 88.652; 
% g=600;
% coords(1:g,1) = coords(1:g,1).*resx;
% coords(1:g,2) = coords(1:g,2).*resy;
% 
% meanpos = mean(coords);
% posx = (coords(1:g,1)-meanpos(1));
% posy = (coords(1:g,2)-meanpos(2));
% varx = var(posx);
% vary = var(posy); 
% kb = 1.38*10^(-5); %in pN*um
% T = 310;
% kx = ((kb*T)/(varx)).*10 %pN/um
% ky = ((kb*T)/(vary)).*10 %pN/um
% 
% KX1=[KX1;kx];
% KY1=[KY1;ky];
% coords2(1:g,1) = coords2(1:g,1).*resx;
% coords2(1:g,2) = coords2(1:g,2).*resy;
% 
% meanpos2 = mean(coords2);
% posx2 = (coords2(1:g,1)-meanpos2(1));
% posy2 = (coords2(1:g,2)-meanpos2(2));
% varx2 = var(posx2);
% vary2 = var(posy2); 
% kb = 1.38*10^(-5); %in pN*um
% T = 310;
% kx2 = ((kb*T)/(varx2)).*10 %pN/um
% ky2 = ((kb*T)/(vary2)).*10 %pN/um
% KX2=[KX2;kx2];
% KY2=[KY2;ky2];
% % 
end
