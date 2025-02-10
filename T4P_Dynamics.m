    %% With this program, you can alyse pilus dynamics - elongation and retraction velocities, from .nd2 data of fluorescently labelled T4P
    close all;
    clc;
    clear;

    %% initiate pathways of folders of your data and the pathway where you want to save your analysed data
    allpath = {'X:\Sebastian\18.03.2022\T136C'}; % Path of your folder with videos

    allsavepath = {'X:\Sebastian\AUSWERTUNG_sebastian\T136C\18.03.2022'};

    %%
    % reference histogram of pixels, to ensure for the right image quality
    datafiles = dir(fullfile('histio.mat'));
    % read in your data
    for ii = 1:length(datafiles)
        tmp = load(fullfile(datafiles(ii).folder, datafiles(ii).name));

        datastr = sprintf('Contour_%u', ii);  % Generate data string
        data.(datastr) = tmp;
    end
    % bring reference histogram in right format
    Bildi=cell2mat(struct2cell(data));
    helpimage=Bildi.Ibild;



    %% analysis of the video



    for l = 1:length(allpath)
        path = allpath{l};
        savepath = allsavepath{l}
        C = strsplit(path,'\');
        Folder = C{end};
        files = dir(fullfile(path,'*.nd2'));
        for d=1:size(files,1)
            close all;
            arrtmp = split(files(d).name,'.');
            arr = arrtmp{1,1};
            arr2 = [Folder,arr];
            if ~exist([savepath,'\',arr2],'dir')
                mkdir([savepath,'\',arr2]);
                xyscale = 0.0792354;
                meta = imreadND2_meta([path '\' arr '.nd2']);
                [img] = imreadND2A( meta,1,1);
                img = squeeze(img);
                allmin = min(img(:));
                allminhelp=allmin;
                allmax = max(img(:));
                allmaxhelp=allmax;
                img(1,1,1) = allmin;
                img=uint8(fix(((1500.0+1)*single(img-allmin)-1)/single(allmax-allmin)));
                %% Choose your cell you want to analyse from overview region of interest (whole field of view)
                f1=figure('Name','Where is the cell?'); imagesc(imhistmatch(img(:,:,1),helpimage)); colormap('winter')
                coord=round(ginput(1)); close(f1)

                %% Set up a region of interest with the right size of your bacterium you want to analyse
                Roisize=40;
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



                % crop image
                imgraw=img(borderx1:borderx2,bordery1:bordery2,:);
                length1=size(imgraw,1);
                length2=size(imgraw,2);
                img=imgraw;


                %%Happy with image quality ?


                sichtbar={'12' '50'};

                while ~isempty(sichtbar)

                    %img wird geändert
                    prompt = {'Happy with video quality?  Range: [9-19]','Happy with contrast?'};
                    dlgtitle = 'Frequency cut-off';
                    definput = sichtbar;
                    dims = [1 60];
                    opts.Interpreter = 'tex';
                    sichtbar = inputdlg(prompt,dlgtitle,dims,definput,opts);

                    if ~isempty(sichtbar)

                        % Assign Cut-off Frequency, via loop
                        D0 =str2num(sichtbar{1,1});
                        cont=str2num(sichtbar{2,1});
                        Movie = zeros(size(img,1),size(img,2),size(img,3),'uint8');
                        %%
                        %This is a new filter which was implemented after publication of 'External Stresses Affect Gonococcal T4P Dynamics' by S Kraus-Römer and I Wieler et al. 2022,Front. Microbiol. 13:839711.
                        %It is a butterworthlowpass filter, allows for dropping high frquency components like noise and by this smoothens the image and enhancing details of the image.
                        for w = 1:size(img,3)


                            input_image=imhistmatch(img(:,:,w),helpimage+cont);

                            [M, N] = size(input_image);

                            % Getting Fourier Transform of the input_image
                            % using MATLAB library function fft2 (2D fast fourier transform)
                            FT_img = fft2(double(input_image));

                            % Assign the order value
                            n = 2; % one can change this value accordingly

                            % Designing filter
                            u = 0:(M-1);
                            v = 0:(N-1);
                            idx = find(u > M/2);
                            u(idx) = u(idx) - M;
                            idy = find(v > N/2);
                            v(idy) = v(idy) - N;

                            % MATLAB library function meshgrid(v, u) returns
                            % 2D grid which contains the coordinates of vectors
                            % v and u. Matrix V with each row is a copy of v
                            % and matrix U with each column is a copy of u
                            [V, U] = meshgrid(v, u);

                            % Calculating Euclidean Distance
                            D = sqrt(U.^2 + V.^2);

                            % determining the filtering mask
                            H = 1./(1 + (D./D0).^(2*n));

                            % Convolution between the Fourier Transformed
                            % image and the mask
                            G = H.*FT_img;

                            % Getting the resultant image by Inverse Fourier Transform
                            % of the convoluted image using MATLAB library function
                            % ifft2 (2D inverse fast fourier transform)
                            output_image = real(ifft2(double(G)));
                            % output_image=uint8(output_image);
                            %
                            %
                            %
                            %  output_image = imadjust(output_image,[0.1 0.5]);
                            Movie(:,:,w)=output_image;

                        end
                        % Display movie of filtered images
                        handle=implay(Movie);
                        mycontrols=handle.DataSource.Controls;
                        play(mycontrols);
                        handle.Parent.Position = [300 300 400 250];
                        set(0,'showHiddenHandles','on');
                        fig_handle = gcf ;
                        fig_handle.findobj % to view all the linked objects with the vision.VideoPlayer
                        ftw = fig_handle.findobj('TooltipString', 'Maintain fit to window');   % this will search the object in the figure which has the respective 'TooltipString' parameter.
                        ftw.ClickedCallback();
                    end
                end

                close(handle)



                % Display movie of final filtered images
                handle=implay(Movie);
                mycontrols=handle.DataSource.Controls;

                handle.Parent.Position = [300 300 400 250];
                set(0,'showHiddenHandles','on');
                fig_handle = gcf ;
                fig_handle.findobj % to view all the linked objects with the vision.VideoPlayer
                ftw = fig_handle.findobj('TooltipString', 'Maintain fit to window');   % this will search the object in the figure which has the respective 'TooltipString' parameter.
                ftw.ClickedCallback();

                %loop if image quality is too bad, but you don't want to stop
                %the loop
                answer = questdlg('Would you like to skip the video?', ...
                    'Skipping',...
                    'yes','no','no');
                % Handle response
                switch answer
                    case 'yes'
                        dessert = 1;
                    case 'no'
                        dessert = 2;
                end

                if dessert==2

                    Moviecutted=Movie(:,:,1:380);
                    Moviecuttedd=imgaussfilt(Moviecutted,5);

                    %% We need the COM for the analysis, here the image will be binarized and center of mass will be determined, if it is not clear where COM is, you can also choose on your own
                    Moviebinarize=imbinarize(Moviecuttedd, 'global'); % binarize image
                    binaryimg=Moviebinarize(:,:,200);

                    CoM=regionprops(Moviebinarize(:,:,200),'centroid'); % Determine center of mass
                    s=CoM;

                    centerofmass = questdlg('Determine center of mass?', ...
                        'Center of mass',...
                        'COM','ME','ME');
                    % Handle response
                    switch centerofmass
                        case 'COM'
                            CoM=regionprops(Moviebinarize(:,:,200),'centroid');
                            s=CoM;
                            cent=s(1).Centroid;
                        case 'ME'
                            f3=figure('Name','Where is the center of the cell?'); imagesc(Movie(:,:,200)); colormap('winter')
                            s=round(ginput(1)); close(f3)
                            cent=s;
                    end

                    %% choose the pili positions which you want to analyse by hand in movie, via drawline selection
                    hold on
                    Steigungenangle=[];
                    answer3=1;
                    count=1;
                    while ~isempty(answer3)
                        count=count+1;
                        prompt = {['Is it the next-to-last Pili you want to analyze??? Choose Cancel!']};
                        dlgtitle = 'Choose LineROi';
                        definput = {'1'};
                        dims = [1 60];
                        opts.Interpreter = 'tex';
                        answer3 = inputdlg(prompt,dlgtitle,dims,definput,opts);

                        %select pili
                        h2 = drawline('SelectedColor','yellow');
                        m2=customWait(h2);
                        hold on;

                        stei=(m2(2,2)-m2(1,2))/(m2(2,1)-m2(1,1));

                        % determine angle, how to rotate the video to analyse
                        % pilus
                        if m2(2,1)<cent(1)
                            if m2(2,2)<cent(2)
                                alpha=180-atand(stei);
                            else
                                alpha=-atand(stei)+180;
                            end
                        else
                            if m2(2,2)<cent(2)
                                alpha=-atand(stei);
                            else
                                alpha=360-atand(stei);
                            end
                        end
                        Steigungenangle=[Steigungenangle; alpha count stei];
                    end



                    %% Choose frames of focus
                    prompt3 = {'Enter startframe:','Enter endframe:','Enter framerate:'};
                    dlgtitle3 = 'Input';
                    dims3 = [1 100];
                    definput3 = {'1','300','19.5'};
                    answer3 = inputdlg(prompt3,dlgtitle3,dims3,definput3)
                    Input3=cellfun(@(x)str2double(x), answer3);
                    frami=Input3(3);
                    framerate=frami;
                    Moviecutted=Movie(:,:,Input3(1):Input3(2));

                    % Arrays which you fill during analysis
                    Langen=[]; % to analyse length of pili
                    Pausingframes=[]; % to analyse pausing of pili
                    Elodurations=[]; % to analyse the duration of elongation
                    Retrdurations=[]; % to analyse the duration of retraction
                    Elovelocities=[]; % to analyse elongation velocities
                    Retrvelocities=[]; % to analyse retraction velocities
                    for a=1:count-1
                        %% rotate video with chosen angle
                        imrot = imrotate(Moviecutted,360-Steigungenangle(a,1),'nearest','crop');

                        %%chose Roi around Pili and cut video
                        handle2 = implay(imrot);
                        handle2.Parent.Position = [100 100 700 550];
                        set(0,'showHiddenHandles','on')
                        fig_handle = gcf ;
                        fig_handle.findobj % to view all the linked objects with the vision.VideoPlayer
                        ftw = fig_handle.findobj ('TooltipString', 'Maintain fit to window');   % this will search the object in the figure which has the respective 'TooltipString' parameter.
                        ftw.ClickedCallback()

                        re = drawrectangle('Position',[cent(1) cent(2) 30 30]);
                        posroi = customWait(re);
                        rectroi=[round(posroi(1)) round(posroi(2)) round(posroi(3)) round(posroi(4))]
                        PiliRoi=imcrop(imrot(:,:,1),rectroi);

                        f3=figure;
                        figure(f3)
                        handle2 = implay(PiliRoi);
                        handle2.Parent.Position = [100 100 700 550];
                        set(0,'showHiddenHandles','on')
                        fig_handle = gcf ;
                        fig_handle.findobj % to view all the linked objects with the vision.VideoPlayer
                        ftw = fig_handle.findobj ('TooltipString', 'Maintain fit to window');   % this will search the object in the figure which has the respective 'TooltipString' parameter.
                        ftw.ClickedCallback()  % ex
                        coord3=round(ginput(1)); close(f3)
                        hold off;
                        % determine intensity profile along the pilus
                        kymo=[];
                        for p =1:size(Moviecutted,3)
                            PiliRoi=imcrop(imrot(:,:,p),rectroi);
                            avarageoverpili=mean(PiliRoi(coord3(2)-2:coord3(2)+2,:)); % intensity profile of 4 Pixels will be averages to ensure that pilus is visible in kymograph later
                            kymo=[kymo; avarageoverpili];
                        end
                        %from the average intensity profiles along the pili,
                        %create a kymograph, the kymograph is also filtered
                        %via wiener2 and gaussian filtering to smooth it. This
                        %helps to determine the edge afterwards, cause we are
                        %interested in the 'contour' of the kymograph which
                        %displays the pilus tip track
                        kymoh=imgaussfilt(wiener2(kymo',[2 2]),3);
                        ia= mat2gray(kymoh(coord3(1):end,:));
                        mina      = max(ia(:));
                        maxa      = min(ia(:));
                        grayUINT8 = uint8((ia - mina) * 255 / (maxa - mina));
                        % figure to check the kymograph
                        figure();
                        kymograph=imshow(grayUINT8);
                        hold off
                        %% evaluate the kymograph by determine the right contrast to detect the edge
                        BW22=grayUINT8;
                        answer1 = 1;
                        imgcopy=BW22;
                        while ~isempty(answer1)
                            %img wird geändert
                            prompt = {'Enter a value for lower limit of contrast (0.01 bis 1): (are you happy ? Click the cancel button!) ', 'Enter value for upper limit (0.01-1), must be greater than lower bound'};
                            dlgtitle = 'Contrast kymograph';
                            definput = {'0.8','0.88'};
                            dims = [1 60];
                            opts.Interpreter = 'tex';
                            answer1 = inputdlg(prompt,dlgtitle,dims,definput,opts);
                            if ~isempty(answer1)
                                bottom=str2num(answer1{1,1});
                                upper=str2num(answer1{2,1});
                                fkymo=figure('Name','Is kymo contrast ok?'); imshow(imadjust(grayUINT8,[bottom upper]));
                            end
                        end
                        close(fkymo);
                        grayUINT8=imadjust(grayUINT8,[bottom upper]);
                        BW22=imgcopy;
                        clear imgcopy;
                        BW2=grayUINT8;
                        answer = 1;
                        imgcopy=BW2;
                        while ~isempty(answer)
                            %img wird geändert
                            prompt = {['Enter a value for sensitivity (0.8 bis 0.01): (are you happy ? Click the cancel button!) ']};
                            dlgtitle = 'Sensitivity Edge detction';
                            definput = {'0.05'};
                            dims = [1 60];
                            opts.Interpreter = 'tex';
                            answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
                            if ~isempty(answer)
                                sensitivity=str2num(answer{1,1});
                                fkymo2=figure('Name','Is kymo ok?'); imshow(edge(wiener2(grayUINT8, [3 3]),'canny',sensitivity,1)); colormap('gray')
                            end
                        end
                        close(fkymo2);
                        BW2=imgcopy;
                        clear imgcopy;
                        BW2 = edge(grayUINT8,'canny',sensitivity,1); % if you are happy with contrast, you can determine the edge via the edge function of matlab
                        BW2=double(BW2);
                        BWd=rmoutliers(BW2);
                        %figure to test the edge
                        figure();
                        ook=imshow(BW2)
                        brush('on')

                        % clear the edge detection and make sure that there is
                        % only one positon of the tip for every time frame
                        for lr=1:size(BW2,2)
                            count=0;
                            for uo=1:size(BW2,1)
                                ou=size(BW2,1)+1-uo;
                                if BW2(ou,lr)==1 && count==0
                                    count=count+1;
                                else
                                    BW2(ou,lr)=0;
                                end
                            end
                        end
                        % test kymograph after clearance
                        figure();
                        ook=imshow(BW2);
                        brush('on')

                        hold off;
                        while size(findobj(ook))>0
                            pause
                        end

                        % get the indices of the kymograph
                        Indizes=[];
                        for de=1:size(BW2,2)
                            for ed=1:size(BW2,1)
                                if BW2(ed,de)==1
                                    Indizes=[Indizes; ed de];
                                end
                            end
                        end

                        %determine surface of bacterium by the minimum of the
                        %pilus, you can directly choose it in the image
                        f6=figure('Name','Where is the background?'); imagesc(BW2); colormap('winter')
                        coord3=round(ginput(1)); close(f6)
                        cellenvelope=coord3(2);

                        % determine edge track substracted by cellenvelope
                        Indi=[]
                        for in=1:length(Indizes(:,1))
                            if Indizes(in,1)>cellenvelope
                                Indi=[Indi;Indizes(in,:)];
                            end
                        end


                        % smooth the edge and check how it looks like
                        figure();
                        peeeek=plot(Indi(:,2),smooth(Indi(:,1)-1-cellenvelope));
                        hold off
                        YDATA=smooth(smooth(Indi(:,1)-1-cellenvelope));
                        XDATA=Indi(:,2);
                        % figure of finale kymograph needs to be closed to
                        % continue with analysis
                        while size(findobj(peeeek))>0
                            pause
                        end

                        %% Now, we analyse retraction velocities, elongation velocities, lengths
                        Masstab= 0.0792354;
                        lengths=findpeaks(Indi(:,1)-cellenvelope,'MinPeakProminence',4)*Masstab; % length are determined via the peaks in the kymopgraph tracks usinf the findpeaks function
                        save('info.mat','framerate','allsavepath','arr2');
                        save('track.mat','XDATA','YDATA');
                        hFig=Dynamicsbearb180222; % GUI to analyse dynamics of T4P

                        set(hFig,'Tag', 'tformChoiceGui','HandleVisibility','on');
                        hGUi=findobj('Tag','formChoiceGui');
                        waitfor(hGUi);
                        while size(findobj(hFig))>0
                            pause
                        end

                        %% Overview
                        % in .mat data of GUI Dynamicsbearb180222 velocities
                        % and lengths will be saved for every single pilus.

                    end
                end


                %% Load data from GUI and put it in array, then save for all analysed pili
                %
                % filename = [[savepath,'\',arr2,'\'],'Pilidynamicsforonebact.mat'];
                % save(filename, 'lengths','Pausingframes','Elodurations','Retrdurations','Elovelocities','Retrvelocities' );
                % close all;


            end
        end


    end










