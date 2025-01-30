% Modified version of bpass. This function hpass attempts to only use the
% high frequency fourier filter used in bpass

% bpass_draft is aiming to let the user choose different band pass filters.
% Optionally low pass and high pass should be able to be chosen from
% different methods
% methods at the moment:
    % 1: Top Head Function Filter in Fourier Space
    % 2: Convolution with 2D gaussian with st.deviation ≃ l / sqrt(2);
    % 3: Convolution with Box Car Average Kernel of circular shape
    % 4: Convolution with Box Car Average Kernel of square shape

function [arr_msk_low_freq,arr_msk_high_freq] = hpass_enno(arr,lobject,method)
if nargin<2
    method=1;
end
res = double(arr);

arr_msk_low_freq = zeros(size(arr));

% circular cut out of low frequencies in fourier space
if method==1

% Object Filters / Low Frequency Filter
    mx = size(arr,2);
    my = size(arr,1);
   
    [dstx2 dsty2] = create_dstx2y2(mx,my);
   
    kx = mx / lobject;
    ky = my / lobject;
    msk_low_freq = ones(my,mx);
    ind = (dstx2./kx^2 + dsty2./ky^2 < 1);
    msk_low_freq(ind) = 0;
    msk_low_freq(floor(my/2) +1,floor(mx/2) +1) = 1;
    
    ft = fft2(res);
    tmp = fftshift(ft);
    tmp3 = tmp .* msk_low_freq;
    arr_msk_low_freq = ifft2(ifftshift(tmp3));
    
% gaussian cut out of low frequencies
elseif method==2
    b = double(lobject);
    
    mx = size(arr,2);
    my = size(arr,1);
    gxy = zeros(my,mx);
    gxy(floor(my/2)+1,floor(mx/2)+1)=1;
    
    kx = mx / double(lobject);
    ky = my / double(lobject);
    
    b = kx/2;
    sm = 0:mx-1;
    r = (sm - floor(mx/2))/(2*b);
    gx = exp(-r.^2)/(2*b*sqrt(pi));
    
    b = ky/2;
    sm = 0:my-1;
    r = (sm - floor(my/2))/(2*b);
    gy = (exp(-r.^2)/(2*b*sqrt(pi)))';
    
    gxy = conv2(gxy,gx,'same');
    gxy = conv2(gxy,gy,'same');
    msk_low_freq = -(gxy./max(gxy(:)))+1;
    msk_low_freq(floor(my/2) +1,floor(mx/2) +1) = 1;
    
    ft = fft2(res);
    tmp = fftshift(ft);
    tmp3 = tmp .* msk_low_freq;
    arr_msk_low_freq = ifft2(ifftshift(tmp3));
    
    msk_high_freq=gxy./max(gxy(:));
    tmp3 = tmp .* msk_high_freq;
    arr_msk_high_freq = ifft2(ifftshift(tmp3));
    
end

  
   
function [dstx2 dsty2] = create_dstx2y2(mx,my)

    centx = floor(mx/2) + 1;
    centy = floor(my/2) +1;
    
    dstx2 = zeros(my,mx);
    for i=1:mx
        dstx2(:,i) = (i-centx).^2;
    end
    dsty2 = zeros(my,mx);
    for i=1:my
        dsty2(i,:) = (i-centy).^2;
    end

