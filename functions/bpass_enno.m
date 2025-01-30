% bpass_draft is aiming to let the user choose different band pass filters.
% Optionally low pass and high pass should be able to be chosen from
% different methods
% methods at the moment:
    % 1: Top Head Function Filter in Fourier Space
    % 2: Convolution with 2D gaussian with st.deviation ≃ l / sqrt(2);
    % 3: Convolution with Box Car Average Kernel of circular shape
    % 4: Convolution with Box Car Average Kernel of square shape
    % 5: masking with gaussian in fourier space

function [res arr_res arr_g w arr_msk_low_freq arr_only_high_freq] = bpass_enno(arr,lnoise,lobject,mthd_noise,mthd_object)

arr = double(arr);

arr_res = zeros(size(arr));
arr_g = zeros(size(arr));
arr_msk_low_freq = zeros(size(arr));
arr_only_high_freq = zeros(size(arr));

if nargin < 4
    mthd_noise = 2;
end
if nargin < 5
    mthd_object = 4;
end

% Noise Filters / High frequency filters
if mthd_noise ==1
    mx = size(arr,2);
    my = size(arr,1);
    [dstx2 dsty2] = create_dstx2y2(mx,my);
    kx = mx / lnoise;
    ky = my / lnoise;
    msk_high_freq = zeros(my,mx);
    ind = (dstx2./kx^2 + dsty2./ky^2 <= 1);
    msk_high_freq(ind) = 1;
    
    res = arr;
%     res = double(DampEdge(res,0.1,1));
    ft = fft2(res);
    tmp = fftshift(ft);
    tmp2 = ((tmp) .* msk_high_freq);
    arr_g = ifft2(ifftshift(tmp2));
    
    msk_only_high_freq = (msk_high_freq -1).* (-1);
    msk_only_high_freq(floor(my/2) +1,floor(mx/2) +1) = 1;
    tmp5 = tmp .* msk_only_high_freq;
    arr_only_high_freq = ifft2(ifftshift(tmp5));
    
    w1 = 0;
end

if mthd_noise == 5
    mx = size(arr,2);
    my = size(arr,1);
    gxy = zeros(my,mx);
    gxy(floor(my/2)+1,floor(mx/2)+1)=1;
    
    kx = mx / double(lnoise);
    ky = my / double(lnoise);
    
    b = kx;
    sm = 0:mx-1;
    r = (sm - floor(mx/2))/(2*b);
%     gx = exp(-r.^2)/(2*b*sqrt(pi));
%     gx = gx./sum(gx(:));
    gx = exp(-r.^2);
    
    b = ky;
    sm = 0:my-1;
    r = (sm - floor(my/2))/(2*b);
%     gy = (exp(-r.^2)/(2*b*sqrt(pi)))';
%     gy = gy./sum(gy(:));
    gy = (exp(-r.^2))';
    
    gxy = conv2(gxy,gx,'same');
    gxy = conv2(gxy,gy,'same');
    msk_high_freq = gxy./max(gxy(:));
%     msk_low_freq(floor(my/2) +1,floor(mx/2) +1) = 1;
    
    res = arr;
%     res = double(DampEdge(res,0.1,1));
    ft = fft2(res);
    tmp = fftshift(ft);
    tmp2 = ((tmp) .* msk_high_freq);
    arr_g = ifft2(ifftshift(tmp2));
    
%     msk_only_high_freq = (msk_high_freq -1).* (-1);
    msk_only_high_freq = -(msk_high_freq./max(msk_high_freq(:)))+1;
    msk_only_high_freq(floor(my/2) +1,floor(mx/2) +1) = 1;
    tmp5 = tmp .* msk_only_high_freq;
    arr_only_high_freq = ifft2(ifftshift(tmp5));
    
    w1 = 0;
end

if mthd_noise ~= 1 && mthd_noise ~=5
    
    b = double(lnoise);
    if mthd_noise == 2
        w1 = round(4*b);
        if w1 < 2
            w1 = 2;
        end
    else
        w1 = ceil(b);
    end
    N = 2*w1 + 1;
    res = arr;
    if mthd_noise == 2
        % Gaussian Convolution kernel
        sm = 0:N-1;
        r = (sm - w1)/(2 * b);
        gx = exp( -r.^2) / (2 * b * sqrt(pi));
        %% 2013-12-01 EO: normalise by the sum?
        gx = gx./sum(gx(:));
        gy = gx';
        g = conv2(res,gx,'same');
        g = conv2(g,gy,'same');
    
    elseif mthd_noise == 3
        x = 0:N-1;
        cent = (N-1) / 2;
        x2= (x-cent).^2;
        dst=zeros(N,N);
        for i=1:N
            dst(i,:)=sqrt((i-1-cent)^2+x2);
        end
        ind = dst<=2*b;
        gxy = zeros(N,N);
        gxy(ind) = 1;
        gxy = gxy ./ sum(gxy(:));
        g = conv2(res,gxy,'same');
    elseif mthd_noise == 4
        gx = zeros(1,N) +1/N;
        gy = gx';
        g = conv2(res,gx,'same');
        g = conv2(g,gy,'same');
    end
%     arr_g = zeros(size(arr));
%     arr_g((w+1):end-w,(w+1):end-w) = g; 
    arr_g = g;
    
end

% Object Filters / Low Frequency Filter
if mthd_object ==1
    mx = size(arr,2);
    my = size(arr,1);
    if mthd_noise ~= 1
        [dstx2 dsty2] = create_dstx2y2(mx,my);
    end
    kx = mx / lobject;
    ky = my / lobject;
    msk_low_freq = ones(my,mx);
    ind = (dstx2./kx^2 + dsty2./ky^2 < 1);
    msk_low_freq(ind) = 0;
    % for arr_res give out the image corresponding to low frequencies
    msk_above_low_freq = (msk_low_freq -1).* (-1);
    msk_low_freq(floor(my/2) +1,floor(mx/2) +1) = 1;
    if mthd_noise ~= 1
        res = arr;
%         res = double(DampEdge(res,0.1,1));
        ft = fft2(res);
        tmp = fftshift(ft);
    end
    tmp3 = tmp .* msk_low_freq;
    
    tmp4 = tmp .* msk_above_low_freq;
    arr_msk_low_freq = ifft2(ifftshift(tmp3));
    arr_res=ifft2(ifftshift(tmp4));
    w2 = 0;
end

if mthd_object ==5

    [arr_msk_low_freq,arr_res] = hpass_enno(arr,lobject,2);
    
    w2 = 0;
end

if mthd_object ~=1 && mthd_object ~=5
    
    b = double(lobject);
    if mthd_object == 2
        w2 = round(4*b);
        if w2 < 2
            w2 = 2;
        end
    else
        w2 = ceil(b);
    end
    
    N = 2*w2 + 1;
    
    res = arr;
    
    if mthd_object == 4
        %Boxcar average kernel: background
        bx = zeros(1,N)  + 1/N;
        by = bx';
        arr_res = conv2(res,bx,'same');
        arr_res = conv2(arr_res,by,'same');
    end
    if mthd_object == 2 
    % Replacement of boxcar kernel with larger gaussian:
        r = (sm - w2)/(2 * b);
        bx = exp( -r.^2) / (2 * b * sqrt(pi));
        by = bx';
        arr_res = conv2(res,bx,'same');
        arr_res = conv2(arr_res,by,'same');
        
    end
    
    if mthd_object == 3
        x = 0:N-1;
        cent = (N-1) / 2;
        x2= (x-cent).^2;
        dst=zeros(N,N);
        for i=1:N
            dst(i,:)=sqrt((i-1-cent)^2+x2);
        end
        ind = dst<=cent;
        bxy = zeros(N,N);
        bxy(ind) = 1;
        bxy = bxy ./ sum(bxy(:));
        arr_res = conv2(res,bxy,'same');
    end
    
%     arr_res = zeros(size(arr));
      
%     arr_res((w+1):end-w,(w+1):end-w) = g2;
    
end

if w1 > w2
    w = w1;
else
    w = w2;
end
msk = zeros(size(arr));
msk((w+1):end-w,(w+1):end-w) = 1;
arr_g = arr_g .* msk;
% arr_res = arr_res .* msk;
if mthd_noise == 1 && mthd_object == 1
    msk_high_low_freq = msk_high_freq .* msk_low_freq;
    res = arr;
%     res = double(DampEdge(res,0.1,1));
    ft = fft2(res);
    tmp = fftshift(ft);
    tmp1 = ((tmp) .* msk_high_low_freq);
    res = ifft2(ifftshift(tmp1));
else
    tmp = arr_res .* msk;
    res = arr_g - tmp;
%     % if negative values should be avoided use this:
%     res = max(arr_g-arr_res,0);
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

