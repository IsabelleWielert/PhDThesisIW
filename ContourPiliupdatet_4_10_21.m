close all; 
clc;
clear all; 
 allpath = {'X:\Sebastian\16.03.2022\T126C'};
%  allpath = {'\\STORM\Images\Sebastian\10.05.2021\test'};


 allsavepath = {'X:\Sebastian\test'};
% allsavepath = {'C:\Users\Isabelle\Desktop\test'};

%% Erstelle Maske 
mapx=710;
mapy=710; 

n=80;
x = 0:mapx;
y = 0:mapy;
[xx yy] = meshgrid(x,y);
u = zeros(size(xx));
u(1:mapx/2,1)=1; 
X1=[u(1:mapx/2,1)];
Y1=zeros(length(X1),1);
Xsrotar=[];
Ystorar=[];

anghelp =-180/pi *(4*pi/(2*n)); % mnur für 90 grad rotation und dann muss ich über die 4 ecken gehen 
count=1;
for j=1:n
i=j;

ang =j*anghelp;  % degrees

Xc =mapx/2;  % X center of rotation
Yc = mapy/2;  % Y center of rotation
% Shift X/Y to the rotation center
Xshift = X1;        
Yshift = Y1;
for i=1:length(X1)
         Xsrot =  round(i*(Xshift(i)*cosd(ang) + Yshift(i)*sind(ang)));
         Ysrot =   round(i*(-Xshift(i)*sind(ang) + Yshift(i)*cosd(ang)));

         Xsrotar=[Xsrotar; Xsrot];
         Ystorar=[Ystorar; Ysrot];
         if Ysrot ==0 || Xsrot==0
         Ysrot=[];
         Xsrot=[];
         end
       u(Xsrot+Xc,Ysrot+Yc)=count;

end

count=count+1;

if j==21
    count=count+1;
elseif j==62
    count=count+1;

end
 
end

% Fülle das U mit den Indizes 

% linksoben:

for i =2:mapx/2+1
    for j =3:mapy/2
        if u(i,j)==0 && u(i,j-1)~= 0
            u(i,j)=u(i,j-1);
        elseif u(i,j)==0 && u(i,j-1)== 0
            u(i,j-2)=61;
        end
   
    end
end

%rechtsoben:


for i =1:mapx/2-1
    for j =1:mapy/2
       
        k=mapy-j+1;
        if u(i,k)==0 && u(i,k+1)~= 0
            u(i,k)=u(i,k+1);
         elseif u(i,k)==0 && u(i,k+1)== 0
            u(i,k-2)=21;
        end
       
    end
end

%linksunten: 
for i =mapx/2:mapx
    for j =3:mapy/2
        if u(i,j)==0 && u(i,j-1)~= 0
            u(i,j)=u(i,j-1);
        elseif u(i,j)==0 && u(i,j-1)== 0
            u(i,j-2)=62;
        end
   
    end
end

%rechtsunten: 



for i =mapx/2:mapx
    for j =1:mapy/2
       
        k=mapy-j+1;
        if u(i,k)==0 && u(i,k+1)~= 0
            u(i,k)=u(i,k+1);
         elseif u(i,k)==0 && u(i,k+1)== 0
            u(i,k-2)=20;
        end
       
    end
end

%delte zeros!

for i =1:mapx
    for j =2:mapy/2
        if u(i,j)==0 && u(i,j-1)~= 0
            u(i,j)=u(i,j-1);
       
        end
   
    end
end



 

%% Bearbeite Bilder und berechne Kontur 



for l = 1:length(allpath)

path = allpath{l};
savepath = allsavepath{l};
C = strsplit(path,'\');
Folder = C{end}; 
files = dir(fullfile(path,'*.nd2'));
for d =1:size(files,1)
close all; 

arrtmp = split(files(d).name,'.');
arr = arrtmp{1,1}; 
arr2 = [Folder,arr];
if ~exist([savepath,'\',arr2],'dir')
mkdir([savepath,'\',arr2])

xyscale = 0.0792354;
meta = imreadND2_meta([path '\' arr '.nd2']);
[img] = imreadND2A( meta,1,1);
   img = squeeze(img);
   allmin = min(img(:));
   allminhelp=allmin;
   allmax = max(img(:)); 
   allmaxhelp=allmax;
   img(1,1,1) = allmin;    
   img=uint8(fix(((255.0+1)*single(img-allmin)-1)/single(allmax-allmin)));
 
f1=figure('Name','Where is the cell?'); imagesc(imadjust(img(:,:,1),[0 0.05])); colormap('winter')
coord=round(ginput(1)); close(f1)

%delete background 
f2=figure('Name','Where is the background?'); imagesc(imadjust(img(:,:,1),[0 0.05])); colormap('winter')
coord2=round(ginput(1)); close(f2)





%% Roi aussuchen und Background abziehen, Kontrast erhöhen
Roisize=50;
if coord(2)-Roisize >0 
    borderx1=coord(2)-Roisize;
else
    borderx1=1; 
end 

if coord(2)+Roisize <size(img,2)
    borderx2=coord(2)+Roisize;
else 
    borderx2=size(img,2); 
end

if coord(1)-Roisize>0
    bordery1=coord(1)-Roisize;
else
    bordery1=1;
end

if coord(1)+Roisize <size(img,2)
    bordery2=coord(1)+Roisize;
else
    bordery2=size(img,2);
end
                



     imgraw=img(borderx1:borderx2,bordery1:bordery2,:);
     length1=size(imgraw,1);
     length2=size(imgraw,2);
     background=img(coord2(2)-floor(length1/2):coord2(2)+floor(length1/2),coord2(1)-floor(length2/2):coord2(1)+floor(length2/2),:);
     img=imgraw-background;
% 
%      figure(); 
%      imshow(img(:,:,1));
%      
%      hold off;

imghelp=imgaussfilt(img(:,:,1));
%  imghelp=imadjust(imghelp,[0 0.8]);
 
%  
% figure(); 
% imshow(imghelp); 
% hold off; 
% imbinarize(I,'adaptive','ForegroundPolarity','dark','Sensitivity',0.4);
for sd=1:size(imghelp,2)
    for ds=1:size(imghelp,1)
        if imghelp(ds,sd)<20
            imghelp(ds,sd)=0; 
        end
    end
end 
binaryimg=imbinarize(imghelp, 'global');

% 
% figure();
% imshow(binaryimg);
% hold off; 
% 


 
Size=9; %control=9; azi=7; cef=13
SE = strel('disk',Size-3,4);
BW2 = imdilate(binaryimg,SE);
SE2 = strel('disk',Size-1,4);
BW3 = imdilate(binaryimg,SE2); 
Mask=BW3-BW2;
% figure();
% imagesc(Mask)
% hold off; 
% %% 

% figure();
% 
% % % imagesc(Mask)
% R(:,:) = uint8(Mask*255);
% G(:,:) =uint8(imghelp);
% B(:,:) = uint8(zeros(size(imghelp)));
% RGB=cat(3,R,G,B);
% imshow(RGB);
% hold off


     
     fltr4img = [1 1 1 1 1; 1 2 2 2 1; 1 2 4 2 1; 1 2 2 2 1; 1 1 1 1 1];
     fltr4img = fltr4img / sum(fltr4img(:));


     Movie = zeros(size(img,1),size(img,2),size(img,3),'uint8');
%      Movie2 = zeros(size(img,1),size(img,2),size(img,3),'uint8');
    for i = 1:size(img,3)
        
    imgunfiltered=wiener2(imadjust(medfilt2(img(:,:,i)),[0 0.02]),[3 3]);% Entscheidend um Kontrast zu beeimnflussen 
    imgfltrd = filter2( fltr4img ,imsharpen(imgunfiltered,'Radius',0.2,'Amount',2) );
    
%     imgunfiltered2=imadjust(img(:,:,i),[0 0.03]);
%     imgfltrd2=filter2( fltr4img , imgunfiltered2 );
%     Movie2(:,:,i)=imgfltrd2;
    Movie(:,:,i)=imgfltrd;
    end
    
    
    
% handle2=implay(Movie2);
% handle2.Parent.Position = [100 100 700 450];
% set(0,'showHiddenHandles','on')
% fig_handle = gcf ;  
% fig_handle.findobj % to view all the linked objects with the vision.VideoPlayer
% ftw = fig_handle.findobj ('TooltipString', 'Maintain fit to window');   % this will search the object in the figure which has the respective 'TooltipString' parameter.
% ftw.ClickedCallback()  
% hold off;    
% 
%  
handle=implay(Movie);
handle.Parent.Position = [100 100 700 550];
set(0,'showHiddenHandles','on')
fig_handle = gcf ;  
fig_handle.findobj % to view all the linked objects with the vision.VideoPlayer
ftw = fig_handle.findobj ('TooltipString', 'Maintain fit to window');   % this will search the object in the figure which has the respective 'TooltipString' parameter.
ftw.ClickedCallback()  

% 
%  while size(findobj(handle))>0
%     pause 
%  end
%   
%  answer = questdlg('Would you like to skip the video?', ...
%      'Skipping',...
%      'yes','no','no');
% % Handle response
% switch answer
%     case 'yes'
%         dessert = 1;
%     case 'no'
%        
%         dessert = 2;
% end
%  
% if dessert==2
% 
%  prompt10 = {'Enter frame without pili'};
%  dlgtitle10 = 'Input';
% dims10 = [1 100];
%  definput10 = {'570'};
%  answer10 = inputdlg(prompt10,dlgtitle10,dims10,definput10)  
%  Input10=cellfun(@(x)str2double(x), answer10);
% 



Moviecutted=Movie(:,:,1:300);
Moviecutted=imgaussfilt(Moviecutted,1);

%%
%%nun bearbeite das Bild um die contour zu kriegen!
Moviebinarize=imbinarize(Moviecutted, 'global');
%% 
Helpim=imbinarize(imadjust(Moviecutted(:,:,200),[0 0.5]),'global');
CoM=regionprops(Helpim,'centroid')
if size(CoM,1)==1

Pili=[];
for t=1:size(Moviecutted,3)
    

% 

s = regionprops(Moviebinarize(:,:,t),'centroid');
% Input7=7; %13 gemessen
% SE = strel('disk',Input7-3,4);
% BW2 = imdilate(Moviebinarize(:,:,Input10),SE);
% SE2 = strel('disk',Input7-1,4);
% BW3 = imdilate(Moviebinarize(:,:,Input10),SE2);
% Mask=BW3-BW2;
%imagesc(Mask)
Arraymask=[];
for h=1:size(Mask,1)
    for g=1:size(Mask,2)
        if Mask(h,g)~=0
            Arraymask=[Arraymask; h g];
        end
    end
end

% 

%% vgl bild mit maps
%1)Moviecutted(:,:,1) original Bild 
%2)Mask, wo ich messen muss 
%3) u Binning Maske, start 


%first step: Verschiebe sector maske: 

try 
mapxs=round(mapx+(Roisize-s.Centroid(2)));
mapys=round(mapy+(Roisize-s.Centroid(1)));

Cutted_u=u(mapxs/2-50:mapxs/2+50,mapys/2-50:mapys/2+50);

% imagesc(Cutted_u)
len_binning=size(Moviecutted(:,:,t),1);
ind_binning=unique(Cutted_u);


    Binnedint=[];
 for a=1:length(ind_binning)
     indi=ind_binning(a);
    help=[];

for g=1:size(Moviecutted(:,:,1),1)
    for f=1:size(Moviecutted(:,:,1),2)

        
         if Mask(g,f)==1 && Cutted_u(g,f)==indi
           
        help=[help; Moviecutted(g,f,t)];
        
        end
        
        
    end
end

    Int_binned=sum(help)/length(help);
  

    Binnedint=[Binnedint; Int_binned] ;
 end
 
 Pili=[Pili Binnedint];


 catch
end
end
% figure(); 
% imagesc(Pili);
% hold off 
% figure(); 
% Pili=255.*Pili./max(Pili);
% imagesc(Pili)
% hold off
% N = 80;                                                         % Number Of Segments
% a = linspace(0, 2*pi, N*10);
% r = 100;
% x = r*cos(a);
% y = r*sin(a);
% hFig=figure;
% 
% plot([zeros(1,N)+50; x(1:10:end)+50], [zeros(1,N)+50; y(1:10:end)+50])
% color = get(hFig,'Color');
% set(gca,'XColor','none','YColor','none','TickDir','out')
% 
% axis tight
% axis equal
% hold off
% F = getframe(gcf);
% [X, Map] = frame2im(F);
% hFig=figure;
% ddd=imshow(X)
% color = get(hFig,'Color');
% set(gca,'XColor','none','YColor','none','TickDir','out')
% dd = ind2rgb(X,Map);
% I = uint8(X,Map);
% I=uint8(I(:,:,1));
% figure();
% imshow(I(290-50:290+50,200-50:200+50,3))
% 
% hFig=figure;
% plot([zeros(1,N)+50; x(1:10:end)+50], [zeros(1,N)+50; y(1:10:end)+50],'color','w','linewidth',1)
% axis tight; 
% hold on
% hf = image(C);
% colormap('gray');
% uistack(hf,'bottom')
% xlim([0 100])
% ylim([0 100])
%  color = get(hFig,'Color');
%  set(gca,'XColor','none','YColor','none','TickDir','out')
%  hold on 
%  rectangle('Position',[80 90 15 4], 'FaceColor', [1 1 1],'EdgeColor',[1 1 1])
% hold off
% 
% C = imfuse(Moviecutted(:,:,100),Mask*255,'falsecolor','Scaling','joint','ColorChannels',[2 1 0]);
% imshow(C)
% 
% figure(); 
% imagesc(Pili);
% hold off 
% %Save
figure();

% % imagesc(Mask)
R(:,:) = uint8(Mask*255);
G(:,:) =uint8(Moviecutted(:,:,3));
B(:,:) =uint8(Cutted_u);
RGB=cat(3,R,G,B);
rgbfigure=imshow(RGB);
hold off
%%

Piliadded=[Pili;Pili(1:5,1:end)];
Piliaddedhelp=Piliadded;
% figure();
% Piliaddim=imagesc(Piliaddedhelp);
% hold off 

%Background abziehen 
  Piliadded=Piliadded-mean(Piliadded(1:end,size(Piliadded,2)));
% figure(); imagesc(Piliadded); hold off;
figure();
imagesc(Piliadded);
hold off 
%bereinigen ? 
for q=1:size(Piliadded,1)
    for p=1:size(Piliadded,2)
        if Piliadded(q,p)<40
            Piliadded(q ,p)=0;
        elseif isnan(Piliadded(q,p))
               Piliadded(q ,p)=0;
        end
    end
end
% figure();
% imagesc(Piliadded);
% hold off 
Piliadded=imgaussfilt(fibermetric(wiener2(Piliadded,[3 3]),3,'ObjectPolarity','bright'),2);
figure();
imagesc(Piliadded);
hold off 
Peaks=[];

for k=1:size(Piliadded,2)
    locs=0.1;
   [pks,locs]=findpeaks(Piliadded(:,k));
   
if locs>0.1
    A=[];
    for p=1:length(locs)
    A=[A;locs(p) k]; 
    end
   Peaks=[Peaks; A];
end

   
end

Pilitofill=zeros(size(Piliadded));



for z=1:length(Peaks)
Pilitofill(Peaks(z,1),Peaks(z,2))=1;
end
% figure();
% imagesc(Pilitofill);
% hold off;

[L,n] = bwlabel(Pilitofill);


Pii=zeros(size(L));

for b=1:n
    [row,col]=find(L==b);
    rowcol=[row col];
    rowstar=rowcol(1,1);
if size(rowcol,1)>10
        for y=rowcol(:,2)
            Pii(rowstar,y)=1;
        end
end
end
% figure();
% h=imagesc(Pii)

Pili=Pii;
Pili=Pili';




figure(); 
handleimgpili=imagesc(Pili);
colormap jet

hold off; 

% 
% C = imfuse(Pii(1:80,1:150).*300,Piliaddedhelp(1:80,1:150),'falsecolor','Scaling','joint','ColorChannels',[2 1 0]);
% imshow(C)
% 



 while size(findobj(handleimgpili))>0
    pause 
 end
  








 
 answer = questdlg('Would you like to save?', ...
     'Skipping',...
     'yes','no','no');
% Handle response
switch answer
    case 'yes'
       savingbutton = 1;
    case 'no'
       
       savingbutton = 2;
end
 




if savingbutton ==1 
    

% 
% filenamefigure='C:\Users\Isabelle\Documents\TestAuswertungPilidynamics\ControlDNA2010\1\Contour.fig';
filenamefigure = [[savepath,'\',arr2,'\'],'rgb.png'];
saveas(rgbfigure, filenamefigure);
filename = [[savepath,'\',arr2,'\'],'Contour.mat'];
save(filename, 'Pili' );
close all; 

else
close all; 
end
end



end
end

end



