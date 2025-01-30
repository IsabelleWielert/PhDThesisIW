clear

%% Parameters
% arr  = 'images\WBmove001.nd2';
% arr  = '3hwt014';
% path = 'D:\Antibiotic project spinning disk\20180928';
% arr  = 'dG4_dPilT_gfp_PI_10xAMP_012';
% path = 'D:\Antibiotic project spinning disk\20181009';
% arr  = 'WT_5ug_ml_cef_3h_002';
% path = 'D:\Antibiotic project spinning disk\20181105';
% arr  = 'wt_8µg_ml_azithromycin_010';
% allpath = {'P:\Antibiotic Project\Control\1h\1h_growth_control', 'P:\Antibiotic Project\Control\1h\20190425', ...
%     'P:\Antibiotic Project\Control\1h\Control 1h', 'P:\Antibiotic Project\Control\4h\20190522', ...
%     'P:\Antibiotic Project\Control\4h\20190523', 'P:\Antibiotic Project\Control\4h\20190528', ... 
%     'P:\Antibiotic Project\Control\6h\20190605', 'P:\Antibiotic Project\Control\6h\20190607', ...
%     'P:\Antibiotic Project\Control\6h\20190612'};
% allsavepath = {'E:\Antibiotic Data\1h Control', 'E:\Antibiotic Data\1h Control', ...
%     'E:\Antibiotic Data\1h Control', 'E:\Antibiotic Data\4h Control', ...
%     'E:\Antibiotic Data\4h Control', 'E:\Antibiotic Data\4h Control', ... 
%     'E:\Antibiotic Data\6h Control', 'E:\Antibiotic Data\6h Control', ...
%     'E:\Antibiotic Data\6h Control'};

% allpath = { 'P:\Antibiotic Project\Control\4hnew\4hControl2', ...
%     'P:\Antibiotic Project\Control\4hnew\20190813', 'P:\Antibiotic Project\Control\4hnew\20190816', ... 
%     'P:\Antibiotic Project\Control\6hnew\6hControl2','P:\Antibiotic Project\Control\6hnew\20190605', 'P:\Antibiotic Project\Control\6hnew\20190607', ...
%     'P:\Antibiotic Project\Control\6hnew\20190612'};
% allsavepath = { 'E:\Antibiotic Data\4h Control cut', ...
%     'E:\Antibiotic Data\4h Control cut', 'E:\Antibiotic Data\4h Control cut', ... 
%     'E:\Antibiotic Data\6h Control cut', 'E:\Antibiotic Data\6h Control cut', ...
%     'E:\Antibiotic Data\6h Control cut','E:\Antibiotic Data\6h Control cut'};

% 
% allpath = { 'P:\Antibiotic Project\Azi\3h\20190207', 'P:\Antibiotic Project\Azi\3h\20190213',...
%     'P:\Antibiotic Project\Azi\3h\20190214', 'P:\Antibiotic Project\Azi\5h\20190702', ... 
%     'P:\Antibiotic Project\Azi\5h\20190703','P:\Antibiotic Project\Azi\5h\azi 6h_alt', ...
%     'P:\Antibiotic Project\Cef\3h Wiederholung\20191107_3hCef', 'P:\Antibiotic Project\Cef\3h Wiederholung\20191108_3hCef', ...
%     'P:\Antibiotic Project\Cef\3h Wiederholung\20191114_3hCef', 'P:\Antibiotic Project\Cef\5h\20190614', ...
%     'P:\Antibiotic Project\Cef\5h\20190625', 'P:\Antibiotic Project\Cef\5h\20190627'};
% 
% allsavepath = { 'E:\Antibiotic Data\3h Azi larger separation2', 'E:\Antibiotic Data\3h Azi larger separation2',...
%     'E:\Antibiotic Data\3h Azi larger separation2', 'E:\Antibiotic Data\5h Azi larger separation2', ... 
%     'E:\Antibiotic Data\5h Azi larger separation2','E:\Antibiotic Data\5h Azi larger separation2', ...
%     'E:\Antibiotic Data\3h Cef larger separation4', 'E:\Antibiotic Data\3h Cef larger separation4', ...
%     'E:\Antibiotic Data\3h Cef larger separation4', 'E:\Antibiotic Data\5h Cef larger separation2', ...
%     'E:\Antibiotic Data\5h Cef larger separation2', 'E:\Antibiotic Data\5h Cef larger separation2'};

% allpath = {'P:\Antibiotic Project\Cef\3h Wiederholung\20191107_3hCef', 'P:\Antibiotic Project\Cef\3h Wiederholung\20191108_3hCef', ...
%     'P:\Antibiotic Project\Cef\3h Wiederholung\20191114_3hCef'};
% allsavepath = {'E:\Antibiotic Data\3h Cef larger separation3', 'E:\Antibiotic Data\3h Cef larger separation3', ...
%     'E:\Antibiotic Data\3h Cef larger separation3'};

% allpath = {'P:\Antibiotic Project\DpilT\5h Cef\20191120_PilT5hCef','P:\Antibiotic Project\DpilT\5h Cef\20191129_PilT5hCef', ...
%     'P:\Antibiotic Project\DpilT\5h Cef\20191130_PilT5hCef'};
% 
% allsavepath = {'E:\Antibiotic Data\PilT 5h CEF large sep','E:\Antibiotic Data\PilT 5h CEF large sep','E:\Antibiotic Data\PilT 5h CEF large sep'};
allpath = {'P:\Antibiotic Project\WB2\3h Cef\20190927_WB2_4hCef', 'P:\Antibiotic Project\WB2\3h Cef\20191023_3hCef', ...
    'P:\Antibiotic Project\WB2\3h Cef\20191030_3hcef'};

allsavepath = {'E:\Antibiotic Data\WB2 3h Cef larger separation', 'E:\Antibiotic Data\WB2 3h Cef larger separation', ...
    'E:\Antibiotic Data\WB2 3h Cef larger separation'};

% allpath = { 'P:\Antibiotic Project\WB2\5h Cef\20191008_WB2_5hCef'};
% 
% allsavepath = { 'E:\Antibiotic Data\WB2 5h Cef'};

% allpath = {'P:\Antibiotic Project\WB2\3h Cef\20191030_3hcef', 'P:\Antibiotic Project\WB2\5h Cef\20191030_5hcef'};
% 
% allsavepath = { 'E:\Antibiotic Data\WB2 3h Cef', 'E:\Antibiotic Data\WB2 5h Cef'};
% 
% allpath = { 'P:\Antibiotic Project\WB2\4h Control\20191022_4hWB2Control', ...
%      'P:\Antibiotic Project\WB2\6h Control\20191022_6hWB2Control'};
% 
% allsavepath = { 'E:\Antibiotic Data\WB2 4h Control',...
%      'E:\Antibiotic Data\WB2 6h Control'};
% allpath = {'P:\Antibiotic Project\DpilT\5h Control\20191108_PilT6h'};
% 
% allsavepath = { 'E:\Antibiotic Data\DpliT 5h Control'};

% allpath = {'P:\Antibiotic Project\pptA\5h Control\20200127'};
% 
% allsavepath = {'E:\Antibiotic Data\DpptA 5h Control'};


for l = 1:length(allpath)

%loop through data
% path = 'P:\Antibiotic Project\Azi\3h\20190214';
% savepath = 'E:\Antibiotic Data\3h Azi';
path = allpath{l};
savepath = allsavepath{l};
C = strsplit(path,'\');
Folder = C{end};
files = dir(fullfile(path,'*.nd2'));
for d = 1:size(files,1)
% for d = 1:1

%     path = 'P:\Antibiotic Project\1h_growth_control';
%     if d < 10
%         arr  = ['20190417_G4_1h_Growth_control_00',num2str(d)];
%     else
%         arr  = ['20190417_G4_1h_Growth_control_0',num2str(d)];
%     end
arrtmp = split(files(d).name,'.');
arr = arrtmp{1,1};
arr2 = [Folder,arr];
mkdir([savepath,'\',arr2])

xyscale = 0.0792354;
zscale = 0.2000000;
% % rsl = [4,5];  %[xy,z]
% % lnoise = [6,6,7];
% % diameter=[51,51,51]; %half
% % sep=[27,27,27];      
% % masksz=[29,29,29]; %set
% % masscut=0;
% % threshold=0;

rsl = [2,5];  %[xy,z]
lnoise = [4,4,4];
% lnoise = [4,4,10];

% diameter=[25,25,25]; %for treatment
diameter=[21,21,21]; %for control
sep=[15,15,15];      %19 always
% sep=[25,25,25];      %19 treatment
% sep=[11,11,11];      %19
masksz=[15,15,15]; %set always
% masksz=[11,11,11]; %set
masscut=0;
threshold=0;
% maxsep = 2;
maxsep = [1, 2, 3, 4, 5];

cthresh = 60; % for control
% cthresh = 40;
% cthresh = 20;

%% Tracking (Position)
% i=1;
% j = 2; % red channel
meta = imreadND2_meta([path '\' arr '.nd2']);
% r=[];
r={};
j = 1;
% for j = 1:1
    for i = 1 : meta.tframes
        % i=1
        
            img = imreadND2A( meta,i,j );
            %     img = squeeze(img(:,:,:,1,2)); % red channel only
            img = squeeze(img);
%             [R,h_COM,rsquare] = meta_col(img);
            
            allmin = min(img(:));
            allmax = max(img(:));
            img(1,1,1) = allmin;
            img=uint8(fix(((255.0+1)*single(img-allmin)-1)/single(allmax-allmin)));
            %         img = resolution(img(:,:,:),rsl);
            rsize = [size(img,1)*rsl(1), size(img,2)*rsl(1), size(img,3)*rsl(2)];
            img = imresize3(img, rsize );
    
        
        %     a = randi([min(min(img(:,:,1))),  mean(mean(img(:,:,1)))],size(img,1)+2*diameter(1),size(img,2)+2*diameter(2),size(img,3)+2*diameter(3));
        %     a(diameter(1)+1:diameter(1)+size(img,1), diameter(2)+1:diameter(2)+size(img,2), diameter(3)+1:diameter(3)+size(img,3)) = img;
        %     a(diameter(1)+1:diameter(1)+size(img,1), diameter(2)+1:diameter(2)+size(img,2), 1:diameter(3)) = repmat(img(:,:,1),[1,1,diameter(3)]);
        %     a(diameter(1)+1:diameter(1)+size(img,1), diameter(2)+1:diameter(2)+size(img,2),diameter(3)+size(img,3)+1:end) = repmat(img(:,:,end),[1,1,diameter(3)]);
        
        %b_tmp = filter9split(a, lnoise, diameter); %bpass für gpu
        b_tmp = bpass3dAW(double(img), lnoise, diameter);
        %     clear a
        %     b = b_tmp(diameter(1)+1:size(img,1)-diameter(1), diameter(2)+1: size(img,2)-diameter(2) , diameter(3)+1:diameter(3)+size(img,3)); %Anton
        b = b_tmp(diameter(1)+1:size(img,1)-diameter(1), diameter(2)+1: size(img,2)-diameter(2) , diameter(3)+1:size(img,3)-diameter(3));
        %     b = b_tmp;
%         if j == 1
%             N_res = size(img,1);
%             
%             imgtmp = imreadND2A( meta,i,2 );
%             
%             imgtmp = squeeze(imgtmp);
%             allmin = min(imgtmp(:));
%             allmax = max(imgtmp(:));
%             imgtmp(1,1,1) = allmin;
%             imgtmp=uint8(fix(((255.0+1)*single(imgtmp-allmin)-1)/single(allmax-allmin)));
%             %         img = resolution(img(:,:,:),rsl);
%             rsizetmp = [size(imgtmp,1)*rsl(1), size(imgtmp,2)*rsl(1), size(imgtmp,3)*rsl(2)];
%             imgtmp = imresize3(imgtmp, rsizetmp );
%             
%             %26.02.2020 x,y beschnitten wie bei b
%             cenimg = single(img(diameter(1)+1:size(img,1)-diameter(1), diameter(2)+1: size(img,2)-diameter(2),diameter(3)+1:size(img,3)-diameter(3)))+single(imgtmp(diameter(1)+1:size(img,1)-diameter(1), diameter(2)+1: size(img,2)-diameter(2),diameter(3)+1:size(img,3)-diameter(3)));
%             clear imgtmp resizetmp
%             allmin = min(cenimg(:));
%             allmax = max(cenimg(:));
%             cenimg(1,1,1) = allmin;
%             cenimg=uint8(fix(((255.0+1)*single(cenimg-allmin)-1)/single(allmax-allmin)));
%             
%             [ centroid, convex ] = calc_colony_centroid( cenimg , cthresh,  savepath, arr2,rsl );
%             
%             clear cenimg
%         end
        
        clear img b_tmp
        
        
        allmin = min(b(:));
        allmax = max(b(:));
        b(1,1,1) = allmin;
        b=fix(((255.0+1)*single(b-allmin)-1)/single(allmax-allmin));
        allmin = b(1,1);
        
        %     savetif(b,[arr '\filtered_image',num2str(i),'.tif'])
        
        
        r_tmp = feature3d( b,masksz, size(b), sep, masscut, threshold);
        %     r_tmp(r_tmp(:,4)<2300,:) = [];
%         r_tmp(r_tmp(:,4)<35000,:) = [];
        if j == 1
            r_tmp(r_tmp(:,4)<15000,:) = [];
        elseif j== 2
            r_tmp(r_tmp(:,4)<15000,:) = [];    
        end
        %     r_tmp(:,3) = r_tmp(:,3) *zscale/rsl(2)/xyscale*rsl(1);
        r_tmp(:,7) = i;
        %     r_tmp(:,8) = j;
        %_tmp(:,8) = size(r,1)+1  : size(r,1) + size(r_tmp,1);

        if j == 2
            [ dist, clust ] = dClust3D( r_tmp, maxsep, rsl, xyscale, zscale);
        end
        %     r = [r; r_tmp];
        r{1,j} = r_tmp;
%             cshow_tc( b*100,r_tmp,path,arr,j);
        cshowRGB_tc( b*100,r_tmp,savepath,arr2,j);
        clear r_tmp
        
    end
% end
% [ r_comb ] = CombinePosChannels( r{1}, r{2}, masksz);
%determine size of clusters from cluster matrix
% [ csize ] = clust_size( clust );

%% Normal
%%{
%calculate distance of dead cells to colony centroid
%careful! x,y axis are switched between r and centroid!
% [ diff_g, normdiff_g ] = dist_to_centr( r, centroid, convex, N_res, 1 );
% [ diff_r, normdiff_r ] = dist_to_centr( r, centroid, convex, N_res, 2 );
% [ diffcont_g, normdiffcont_g ] = dist_to_cont( r, centroid, convex, N_res, 1 );
% [ diffcont_r, normdiffcont_r ] = dist_to_cont( r, centroid, convex, N_res, 2 );
% [ diffcont_all, normdiffcont_all ] = dist_to_cont( r_comb, centroid, convex, N_res );
% % [Ng,edgesg] = histcounts(normdiff_g,'BinLimits',[0,ceil(max(normdiff_r)*10)/10],'BinWidth',0.1);
% % [Nr,edgesr] = histcounts(normdiff_r,'BinLimits',[0,ceil(max(normdiff_r)*10)/10],'BinWidth',0.1);
% 
% [Ng,edgesg] = histcounts(normdiff_g,'BinLimits',[0,1.5],'BinWidth',0.1);
% [Nr,edgesr] = histcounts(normdiff_r,'BinLimits',[0,1.5],'BinWidth',0.1);
% 
% Norm_dead = Nr./Ng(1:length(Nr));
% x = [0.05:0.1:edgesr(end-1)+0.05];
% 
% Norm_dead(isnan(Norm_dead))=0;
% Norm_dead(isinf(Norm_dead))=1;

% 
% figure();plot(x,Norm_dead,'+')
% figure();
% histogram('BinEdges',edgesg,'BinCounts',Ng)
% figure();
% histogram('BinEdges',edgesr,'BinCounts',Nr)
% figure();histogram(csize,'BinMethod','Integers');



savefile = [savepath,'\',arr2,'\data.mat'];
save(savefile,'r','N_res')
%}

%% PilT
%{ 
savefile = [savepath,'\',arr2,'\data.mat'];
save(savefile,'r','r_comb','csize','maxsep','cthresh','centroid','convex','N_res')
%}


% cimg = loadtiff([savepath,'\',arr2,'\centroid_images.tif']);
% cshowRGB_centroid_tc( b,10*uint16(cimg), r,savepath,arr2,1,diameter,rsl);
% delete([savepath,'\',arr2,'\centroid_images.tif'])
end
end
% clear b
% cshow(b*100,r)
% sshow(b,r);





%% Anton Tracking etc
%{
%% Tracking (Trajectories)

zz=r;
% r=zz;

r_tmp = r;

%maximal displacement (in pixels) a feature may undergo between successive frames
% goodenough = uint8(meta.tframes) ;   %minimum length requirement for a trajectory to be retained/ bewahrt                         %how many consequtive frames a feature is allowed to skip

goodenough = 5 ;
memory = 1;
maxdisp = 20;

res = trackmem(r, maxdisp,3 ,goodenough ,memory);

size(res)
size(r)

% res=[];  
% for maxdisp = 30
%     res_tmp = trackmem(r, maxdisp,3 ,goodenough ,memory);
%     r_tmp(res_tmp(:,8),:)=0;
%     r=r_tmp;
%     r( r_tmp(:,8) == 0,: ) = [];
%     res_tmp(:,6) = maxdisp;
%     res = [res; res_tmp];
% end


% res = trackmem(r, 15,3 ,2 ,0);

pixel_bias(res);
plot1_traj( res )
subplot(3,3,7)
plot_traj(res);
title('Trajectories 3D')
subplot(3,3,8)
plot(res(:,6),sqrt(res(:,5)),'.')
xlabel('integrated brightness')
ylabel('radius of gyration')
subplot(3,3,9)
plot(res(:,4),res(:,6),'.')
xlabel('integrated brightness')
ylabel('peak height of the feature')
savefig([ path '\' arr '\' arr '.fig'])
% 



[~, pos_tmp] = sort(res(:,7));
res_sort = res(pos_tmp,:);
traj = histcounts(res_sort(:,9),max(res(:,9)));
[traj_sort, pos_tmp2 ] = sort(traj); %von kurz zu lang

% save([ arr '\' arr '_workspace'])


% res(:,1:3) = res(:,1:3)*xyscale/rsl(1);
for i=1:meta.tframes
    idx = find(res(:,7)==i);
    ovitotxt(res(idx,:) ,[path '\' arr '\particle_location_',num2str(i),'.txt'])
end



%Diffusion(radius)

Diff_r=[];
center = [mean(res(:,1)),mean(res(:,2))];
for i = 1 : max(res(:,9))
    pos_tmp = find(res(:,9)==i);
    msd_tmp = MSD(res(pos_tmp,1:3));
    Diff =  mean(sqrt(msd_tmp(:,1)).*msd_tmp(:,3));
    
    dist = sqrt((center(1) - res(pos_tmp(1),1))^2 + (center(2) - res(pos_tmp(1),2))^2);
    Diff_r = [Diff_r ; [dist, Diff]]; 
end
plot(Diff_r(:,1),Diff_r(:,2),'x')



%% autocorrelation function in X
n = 7; %lag number
c = zeros(n,3);
c(1,1:3) = 1;
for i = 1 : max(res(:,9))
    pos_tmp = find(res(:,9)==i);
    if length(pos_tmp) > 15 
        res_tmp = res(pos_tmp,1:3);
        v = res_tmp(2:end,1:3) - res_tmp(1:end-1,1:3);
        %Normalization factors
        c_x = sum( v(1:end,1) .* v(1:end,1) );
        c_y = sum( v(1:end,2) .* v(1:end,2) );
        c_z = sum( v(1:end,3) .* v(1:end,3) );
        for t = 1:n-1
            c(t+1,1) = sum( v(1:end-t,1) .* v(1+t:end,1) )/c_x;
            c(t+1,2) = sum( v(1:end-t,2) .* v(1+t:end,2) )/c_y;
            c(t+1,3) = sum( v(1:end-t,3) .* v(1+t:end,3) )/c_z;
        end
        
        if max(abs(c(2,:)))>0.5
        
            figure
            subplot(3,2,1)
            plot(0:n-1,c(:,1))
            subplot(3,2,3)
            plot(0:n-1,c(:,2))
            subplot(3,2,5)
            plot(0:n-1,c(:,3))

            subplot(3,2,2)
            plot(v(:,1))
            hold on
            plot([0 ; v(:,1)])
            subplot(3,2,4)
            plot(v(:,2))
            hold on
            plot([0 ; v(:,2)])
            subplot(3,2,6)
            plot(v(:,3))
            hold on
            plot([0 ; v(:,3)])
            
        end
    end
end

%% cross correlation

n = 7; %lag number
c = zeros(n,3);
c(1,1:3) = 1;

for i = 1 : max(res(:,9))
    pos_tmp = find(res(:,9)==i);
    
    res_tmp = res(pos_tmp,:);
    time_point = res_tmp(1,7);
    
    pos_time = find(res(:,9) == time_point);
    
    
        
        
        
        
        
        
        
        
        
        
        v = res_tmp(2:end,1:3) - res_tmp(1:end-1,1:3);
        %Normalization factors
        c_x = sum( v(1:end,1) .* v(1:end,1) );
        c_y = sum( v(1:end,2) .* v(1:end,2) );
        c_z = sum( v(1:end,3) .* v(1:end,3) );
        for t = 1:n-1
            c(t+1,1) = sum( v(1:end-t,1) .* v(1+t:end,1) )/c_x;
            c(t+1,2) = sum( v(1:end-t,2) .* v(1+t:end,2) )/c_y;
            c(t+1,3) = sum( v(1:end-t,3) .* v(1+t:end,3) )/c_z;
        end
        
        figure
        subplot(3,2,1)
        plot(0:n-1,c(:,1))
        subplot(3,2,3)
        plot(0:n-1,c(:,2))
        subplot(3,2,5)
        plot(0:n-1,c(:,3))
        
        subplot(3,2,2)
        plot(v(:,1))
        hold on
        plot([0 ; v(:,1)])
        subplot(3,2,4)
        plot(v(:,2))
        hold on
        plot([0 ; v(:,2)])
        subplot(3,2,6)
        plot(v(:,3))
        hold on
        plot([0 ; v(:,3)])
    end
end





%msd in 6st position
for i = 1 : max(res(:,9))
    pos_tmp = find(res(:,9) == i);
    msd_tmp = MSD(res(pos_tmp,:));
    res(pos_tmp,6) = mean(msd_tmp(:,1));
end
%}
