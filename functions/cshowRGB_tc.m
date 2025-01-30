function [ ] = cshowRGB_tc( img  , mx,path,arr,ch)

mask_radius = [15,15,15];

[nx,ny,nz] = size(img);
bignum = nx*ny*nz;

%build large mask and image with masks
extent = fix(mask_radius*2) + 1	;
extent = extent + mod(extent+1, 2);
rsq = lrsqd3dMB(extent, mask_radius(1,3)/mask_radius(1,1));
mask = rsq < (mask_radius(1))^2;
% cast the mask into a one dimensional form-- imask!
bmask = zeros(nx,ny,extent(3));
bmask(1:extent(1),1:extent(2),:) = mask;
imask = find(bmask > 0) + bignum -(nx*ny*(fix(extent(3)/2)))-(nx*fix((extent(2)/2))) -fix(extent(1)/2);


mx=uint16(mx);
r=sub2ind( size(img), mx(:,1), mx(:,2), mx(:,3) );
temp = img;
for i=1:length(r)
    indx= mod(r(i)+imask-2,bignum)+1; 
    ndiff = find(abs(diff(indx))> nx*ny);
    if isempty(ndiff)
        img(indx) = img(indx) + 4000;
    elseif mx(i,3) < mask_radius(1)
        img(indx(ndiff+1:end)) = img(indx(ndiff+1:end)) + 4000;
    elseif mx(i,3) > mask_radius(1)
        img(indx(1:ndiff)) = img(indx(1:ndiff)) + 4000;
    end
%     img(indx) = img(indx) + 4000;
end

for i=1:size(mx,1)
    x=mx(i,1);
    y=mx(i,2);
    z=mx(i,3);
    img(x,y,z)=65535;
%     img(x+1,y,z)=65535;
%     img(x,y+1,z)=65535;
%     img(x-1,y,z)=65535;
%     img(x,y-1,z)=65535;
%     img(x,y,z+1)=65535;
%     img(x,y,z-1)=65535;

%     img(x+1,y+1,z)=65535;
%     img(x+1,y-1,z)=65535;
%     img(x-1,y+1,z)=65535;
%     img(x-1,y-1,z)=65535;
%     
%     img(x+1,y,z-1)=65535;
%     img(x,y+1,z-1)=65535;
%     img(x-1,y,z-1)=65535;
%     img(x,y-1,z-1)=65535;
% 
%     img(x+1,y,z+1)=65535;
%     img(x,y+1,z+1)=65535;
%     img(x-1,y,z+1)=65535;
%     img(x,y-1,z+1)=65535;
%     
%     img(x+2,y,z)=65535;
%     img(x,y+2,z)=65535;
%     img(uint16(x-2),y,z)=65535;
%     img(x,uint16(y-2),z)=65535;
end
RGB = uint16(cat(4, img, temp, temp));
RGB = permute(RGB,[1,2,4,3]);
% savetif(img,['green_filtered_image_',num2str(1),'.tif'])
% imshow3D(img);
% savetif(img,[path, '\', arr, '\filtered_image_mask_', num2str(ch), '.tif'])

options.color = true;
saveastiff(RGB,[path, '\', arr, '\filtered_image_mask_', num2str(ch), '.tif'],options);
end