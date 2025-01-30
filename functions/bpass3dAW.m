function res=bpass3dAW(image, lnoise, lobject)

%bpass3d is written by Yongxiang Gao and Maria Kilfoil, based on
%the IDL code written by John C. Crocker and David G. Grier.
%the input variable supposed to be contained in varargin is [noclip nopad]
%inputv is to indicate whether there is input for noclip or nopad by
%logical number 1 and 0. 
%Last visit is on June 15, 2005

% ;		Implements a real-space bandpass filter to suppress 
% ;		pixel noise and slow-scale image variations while 
% ;		retaining information of a characteristic size.
% ;		*Works with anisotropic 3d cube data*

%Additionally modified by Anton Welker 4.1.2017

[nx,ny,nf]=size(image);

% do xdirection masks
bb=single(lnoise(1));
w = lobject(1);
N = 2*w + 1;
r = single(-w:w)/(2 * bb);
gx = exp( -r.^2 );
gx = gx /sum(gx);
bx = zeros(1,N,'single') - 1./N;
factor = ( sum(gx.^2) - fix(1/N) );
gy =(gx)';
by =(bx)';


            	
% do z direction masks
bb = single(lnoise(2)); 
w = round( max(lobject(3),2 * bb) );
N = 2*w + 1;
r = single(-w:w)/(2 * bb);
gz = exp( -r.^2 );
gz = gz / sum(gz);
gz = (gz)';
bz = zeros(1,N,'single') - 1/N;
bz = (bz)'; 


b=image;
% do x and y convolutions
for i=1:nf
    image(:,:,i)=conv2(gx,gy,image(:,:,i),'same');
    b(:,:,i)=conv2(bx,by,b(:,:,i),'same');
end

% do z convolution
clear temp tep temp2 tep2 d f 

%%{
temp=permute(image, [3,1,2]);
d=zeros(nf,nx,ny);

for i=1:ny
    d(:,:,i)= conv2(temp(:,:,i),gz,'same');
end

clear temp
image=permute(d,[2,3,1]);
clear d
temp2=permute(b,[3,1,2]);
clear b

temp2=gather(temp2);

f=zeros(nf,nx,ny);
for i=1:ny
   f(:,:,i)= conv2(temp2(:,:,i),bz,'same');
end 
clear temp2
b=permute(f,[2,3,1]);
clear d f temp temp2
%}

%{
temp=permute(image, [3,2,1]);
d=zeros(nf,nx,ny);

for i=1:ny
    d(:,:,i)= conv2(temp(:,:,i),gz,'same');
end

clear temp
image=permute(d,[3,2,1]);
clear d
temp2=permute(b,[3,2,1]);
clear b

temp2=gather(temp2);

f=zeros(nf,nx,ny);
for i=1:ny
   f(:,:,i)= conv2(temp2(:,:,i),bz,'same');
end 
clear temp2
b=permute(f,[3,2,1]);
clear d f temp temp2
%}


image=image+b; % This is the final straw on memory: 12X number of original
clear b temp

res=max(image./factor,0);
% res=max(image,0);
res=gather(res);


% %% adjust intesity
% [~,~,nz] = size(res);
% 
% for i=1:3
%     part_max = max(max(max((res(:,:,i)))));
%     res(:,:,i)=255/part_max*res(:,:,i);
% end
% for i=4:nz-3
%     part_max = max(max(max((res(:,:,i-3:i+3)))));
%     res(:,:,i)=255/part_max*res(:,:,i);
% end

