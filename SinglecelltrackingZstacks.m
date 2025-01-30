%% This program was developed by Tom Cronenberg and Marc Hennes and Anton Welker. Parameters were adjusted by me and it was used to track cells in 2D slices of a z-stacks (Spinning Disk Microscopy)
%% Parameters

% Enter your folder names
allpath = {'Folder of your data .nd2'};
allsavepath = {'Save path names'};


for l = 1:length(allpath)
    % open z-stacks
    path = allpath{l};
    savepath = allsavepath{l};
    C = strsplit(path,'\');
    Folder = C{end};
    files = dir(fullfile(path,'*.nd2'));
    for d = 1:size(files,1)
        
        
        arrtmp = split(files(d).name,'.');
        arr = arrtmp{1,1};
        arr2 = [Folder,arr];
        mkdir([savepath,'\',arr2])
        
        % set xy sclae and z-scale
        xyscale = 0.0792354;
        zscale = 0.2000000;
        
        
        rsl = [2,5];  %[xy,z]
        lnoise = [4,4,4];
        
        diameter=[21,21,21]; %for control
        sep=[19,19,19];      %19 always von peak zu peak durchmesser von einer maske, zellgroesse 1px=
        masksz=[15,15,15]; %set always
        
        masscut=0;
        threshold=0;
        
        maxsep = [1, 2, 3, 4, 5];% threshold for separation
        
        cthresh = 60; % for control
        
        
        %% Tracking (Position)
        
        
        meta = imreadND2_meta([path '\' arr '.nd2']);
        r={};
        j = 1;
        % for j = 1:1
        for i = 1 : meta.tframes
            % i=1
            
            
            img = imreadND2A( meta,100,1 );
            
            img = squeeze(img);
            %
            
            allmin = min(img(:));
            allmax = max(img(:));
            img(1,1,1) = allmin;
            img=uint8(fix(((255.0+1)*single(img-allmin)-1)/single(allmax-allmin)));
            
            rsize = [size(img,1)*rsl(1), size(img,2)*rsl(1), size(img,3)*rsl(2)]; % times rsl to get the right sizes
            img = imresize3(img, rsize );
            
            
            b_tmp = bpass3dAW(double(img), lnoise, diameter);
            b = b_tmp(diameter(1)+1:size(img,1)-diameter(1), diameter(2)+1: size(img,2)-diameter(2) , diameter(3)+1:size(img,3)-diameter(3));
            
            
            clear img b_tmp
            
            
            allmin = min(b(:));
            allmax = max(b(:));
            b(1,1,1) = allmin;
            b=fix(((255.0+1)*single(b-allmin)-1)/single(allmax-allmin));
            allmin = b(1,1);
            
            
            r_tmp = feature3d( b,masksz, size(b), sep, masscut, threshold);
            r_tmp(r_tmp(:,4)<15000,:) = [];
            r_tmp(r_tmp(:,3)>210,:) =[];
            r_tmp(:,7) = i;
            
            
            
            
            r{1,j} = r_tmp;
            %
            cshowRGB_tc( b*100,r_tmp,savepath,arr2,j);
            clear r_tmp
            %
        end
        
        %% Save positions of cells
        
        savefile = [savepath,'\',arr2,'\data.mat'];
        save(savefile,'r')
        
        
    end
end

