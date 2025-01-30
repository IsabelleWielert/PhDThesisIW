    %% Program to analyse contour of single cells, determine pilus production rate
    %% The contour will be sectioned in 80 parts via a mask. Then, the image will be processed to increase the contrast.
    %% In the last step the contoouor is created and the pass through this contour in a certain section will be measured.
    %% From these a time kymograph is created in which the passes of the T4P can be counted.

    % for new data, go through code and test your analysis, adjust
    % thresholds, contrast and filters



    close all;
    clear all;
    clc;
    % Read in your file where you saved the videos .nd2 format
    allpath = {'your folder path'};

    % Create path to save you analysis
    allsavepath = {'your folder path'};

    %% Create circular mask for sections
    mapx=710;
    mapy=710;

    n=80 % number of sections

    x = 0:mapx;
    y = 0:mapy;
    [xx yy] = meshgrid(x,y);
    u = zeros(size(xx));
    u(1:mapx/2,1)=1;
    X1=[u(1:mapx/2,1)];
    Y1=zeros(length(X1),1);
    Xsrotar=[];
    Ystorar=[];

    anghelp =-180/pi *(4*pi/(2*n));
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

    % Fill up your mask matrix with indices

    % left upper quartile

    for i =2:mapx/2+1
        for j =3:mapy/2
            if u(i,j)==0 && u(i,j-1)~= 0
                u(i,j)=u(i,j-1);
            elseif u(i,j)==0 && u(i,j-1)== 0
                u(i,j-2)=61;
            end

        end
    end

    %right top quartile


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

    %left bottom quartile
    for i =mapx/2:mapx
        for j =3:mapy/2
            if u(i,j)==0 && u(i,j-1)~= 0
                u(i,j)=u(i,j-1);
            elseif u(i,j)==0 && u(i,j-1)== 0
                u(i,j-2)=62;
            end

        end
    end

    %rigth bottom quartile

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

    %delete zeros

    for i =1:mapx
        for j =2:mapy/2
            if u(i,j)==0 && u(i,j-1)~= 0
                u(i,j)=u(i,j-1);

            end

        end
    end

    % Figure to test the mask:
    figure();
    imagesc(u);
    hold off

    bins=unique(u); % Test if you have unique bins


    %% Pre-process the images


    % read out the file
    for l = 1:length(allpath)

        path = allpath{l};
        savepath = allsavepath{l};
        C = strsplit(path,'\');
        Folder = C{end};
        files = dir(fullfile(path,'*.nd2'));
        for d =1:size(files,1)
            % bring frames of video in right format

            arrtmp = split(files(d).name,'.');
            arr = arrtmp{1,1};
            arr2 = [Folder,arr];
            if ~exist([savepath,'\',arr2],'dir')
                mkdir([savepath,'\',arr2])

                xyscale = 0.0792354; % Adjust xy scale of images
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
                coord=round(ginput(1)); close(f1) % Choose cell in graphical user interface

                %delete background
                f2=figure('Name','Where is the background?'); imagesc(imadjust(img(:,:,1),[0 0.05])); colormap('winter')
                coord2=round(ginput(1)); close(f2) % choose position of background to substract it from image intensity



                %%  Determine ROI and substract background, increase the contrast to create the mask for the contour

                Roisize=50; % Set Roi Size, if you would like to analyse colony, increase
                % loop to crop image:
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

                % end up with img, which is cropped and background is substracted

                imgraw=img(borderx1:borderx2,bordery1:bordery2,:);
                length1=size(imgraw,1);
                length2=size(imgraw,2);
                background=img(coord2(2)-floor(length1/2):coord2(2)+floor(length1/2),coord2(1)-floor(length2/2):coord2(1)+floor(length2/2),:);
                img=imgraw-background;


                % gaussian filter to smooth the background
                imghelp=imgaussfilt(img(:,:,1));
                % dependent on image quality set an intensity threshold I_threshold to discard pixel
                % intensities
                I_threshold=20;
                for sd=1:size(imghelp,2)
                    for ds=1:size(imghelp,1)
                        if imghelp(ds,sd)<I_threshold
                            imghelp(ds,sd)=0;
                        end
                    end
                end

                %% Create contour: first binarize image,

                %binarize image:
                binaryimg=imbinarize(imghelp, 'global');

                %create mask, by first expand binarized image and then substract to end up
                %with contour:

                Size=8; %control=9; azi=7; cef=13 (Size for expanding, depends on cell size, so make sure that this value matches your needs)
                SE = strel('disk',Size-3,4);
                BW2 = imdilate(binaryimg,SE); % expand
                SE2 = strel('disk',Size-1,4);
                BW3 = imdilate(binaryimg,SE2); % expand
                Mask=BW3-BW2; % Finale mask by substracting the two expansions

                %% Process image further to get the pili visible

                fltr4img = [1 1 1 1 1; 1 2 2 2 1; 1 2 4 2 1; 1 2 2 2 1; 1 1 1 1 1]; % Median Filter
                fltr4img = fltr4img / sum(fltr4img(:));


                Movie = zeros(size(img,1),size(img,2),size(img,3),'uint8'); %now loop over the whole movie to process the frames
                for i = 1:size(img,3)
                    imgunfiltered=wiener2(imadjust(medfilt2(img(:,:,i)),[0 0.1]),[4 4]);% set the contrast, depends on image quality %before 0 0.02 and 3 3
                    imgfltrd = filter2( fltr4img ,imsharpen(imgunfiltered,'Radius',0.2,'Amount',1) ); % sharpend and median filter to increase image quality
                    Movie(:,:,i)=imgfltrd;
                end
                % Play movie and control if you are satisfied with the image settings
                handle=implay(Movie);
                handle.Parent.Position = [100 100 700 550];
                set(0,'showHiddenHandles','on')
                fig_handle = gcf ;
                fig_handle.findobj % to view all the linked objects with the vision.VideoPlayer
                ftw = fig_handle.findobj ('TooltipString', 'Maintain fit to window');   % this will search the object in the figure which has the respective 'TooltipString' parameter.
                ftw.ClickedCallback()

                % Cut the movie if necessary, here 380 frames:
                Moviecutted=Movie(:,:,1:380);
                Moviecutted=imgaussfilt(Moviecutted,1); % apply again a gaussian filter to smooth the images:

                % determine center of mass
                Moviebinarize=imbinarize(Moviecutted, 'global');
                CoM=regionprops(Moviebinarize(:,:,200),'centroid'); % is determined from 200th frame here
                Pili=[];

                if size(CoM,1)==1 % just to make sure, there is no neighbouring cell and the Center of Mass is well defined
                    for t=1:size(Moviecutted,3)

                        s = regionprops(Moviebinarize(:,:,t),'centroid'); % determine for every frame is cell is moving
                        % get get indices from conour in array, just another format:
                        Arraymask=[];
                        for h=1:size(Mask,1)
                            for g=1:size(Mask,2)
                                if Mask(h,g)~=0
                                    Arraymask=[Arraymask; h g];
                                end
                            end
                        end
                        % read out pixel values along contour

                        try %sometimes errors happens here and I dont want the whole folder to stop analysing

                            % Crop the movie
                            mapxs=round(mapx+(Roisize-s.Centroid(2))); % crop the movie
                            mapys=round(mapy+(Roisize-s.Centroid(1))); % crop the movie

                            Cutted_u=u(mapxs/2-50:mapxs/2+50,mapys/2-50:mapys/2+50);
                            % now compare mask and image to get intensity values:
                            len_binning=size(Moviecutted(:,:,t),1);
                            ind_binning=unique(Cutted_u);

                            Binnedint=[];
                            for a=1:length(ind_binning)
                                indi=ind_binning(a);
                                help=[];
                                for g=1:size(Moviecutted(:,:,1),1)
                                    for f=1:size(Moviecutted(:,:,1),2)
                                        if Mask(g,f)==1 && Cutted_u(g,f)==indi
                                            help=[help; Moviecutted(g,f,t)]; % sometimes if contour is not smooth more than 1 Pixel beongs to one section, take the average:
                                        end
                                    end
                                end
                                Int_binned=sum(help)/length(help); % For averaging
                                Binnedint=[Binnedint; Int_binned] ;
                            end

                            Pili=[Pili Binnedint]; % Save pixels values of contour for each frame to build kymograph

                        catch
                        end
                    end

                    % check how Pili plot look like:
                    figure();
                    imagesc(Pili);
                    hold off

                    figure();
                    Pili=255.*Pili./max(max(Pili));% just increase the contrast a bit for the contur
                    imagesc(Pili)
                    hold off


                    figure();
                    imagesc(Pili);
                    hold off


                    %%Save RGB Image to see how contour looks like:
                    figure();

                    R(:,:) = uint8(Mask*255);
                    G(:,:) =uint8(Moviecutted(:,:,3));
                    B(:,:) = uint8(Cutted_u.*255./40);
                    RGB=cat(3,R,G,B);
                    rgbfigure=imshow(RGB);
                    hold off
                    %% create finale kymograph
                    Piliadded=[Pili;Pili(1:5,1:end)]; % Extend Matrix if there pili at the edge of (1 or 80 section)
                    figure();
                    imagesc(Piliadded);
                    hold off

                    %substract once again the backrgound of the kymograph to
                    %make analysis easier
                    Piliadded=Piliadded-mean(Piliadded(1:end,size(Piliadded,2)-10:size(Piliadded,2)),2);figure(); imagesc(Piliadded); hold off;
                    % check if this is satisying, otherwise adjust
                    figure();
                    imagesc(Piliadded);
                    hold off
                    %clean up, if intensities are below threshold: here 40
                    I_p_threshold=40;
                    for q=1:size(Piliadded,1)
                        for p=1:size(Piliadded,2)
                            if Piliadded(q,p)<I_p_threshold;
                                Piliadded(q ,p)=0;
                            end
                        end
                    end
                    figure();
                    imagesc(Piliadded);
                    hold off
                    % apply first a wiener filter, than increase filamentous
                    % strcutures and next smooth the kymograph via a gaussian
                    % filter:
                    Piliadded=imgaussfilt(fibermetric(wiener2(Piliadded,[3 3]),3,'ObjectPolarity','bright'),1);

                    % check if the kymograph is satisyfying?
                    figure();
                    imagesc(Piliadded);
                    hold off


                    Peaks=[];
                    for k=1:size(Piliadded,2)
                        locs=0.1;
                        % findpeaks is a matlab function, here useful to find
                        % the pili along the contour
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
                        Pilitofill(Peaks(z,1),Peaks(z,2))=1; % create binary 'where is pili present' mask
                    end

                    [L,n] = bwlabel(Pilitofill);
                    Pii=zeros(size(L));
                    %loop to clean up matrix of where pili are present
                    for b=1:n
                        [row,col]=find(L==b);
                        rowcol=[row col];
                        rowstar=rowcol(1,1);
                        if size(rowcol,1)>8
                            for y=rowcol(:,2)
                                Pii(rowstar,y)=1;
                            end
                        end
                    end

                    %just transform matrix
                    Pili=Pii;
                    Pili=Pili';



                    %final figure for where are pili present!
                    figure();

                    Red(:,:) = uint8(Pili'*150);
                    Blue(:,:) = uint8(zeros(size(Pili',1),(size(Pili',2))));
                    Green(:,:) =uint8(Piliadded*100);
                    RGB1=cat(3,Red,Green,Blue);
                    handleimgpili=imshow(RGB1,'InitialMagnification',1000);
                    hold off

                    while size(findobj(handleimgpili))>0
                        pause
                    end


                    %% Here GUI to check if you are satisfied with your analysis:

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
                        % Here yo can enter name of rgb image and the file
                        % where you saved the contour
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

