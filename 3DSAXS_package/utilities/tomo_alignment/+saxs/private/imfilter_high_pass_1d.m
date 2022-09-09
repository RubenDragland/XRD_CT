function img = imfilter_high_pass_1d(img, ax, sigma)

    % IMFILTER_HIGH_PASS_1D applies fft filter along AX dimension that
    % removes SIGMA ratio of the low frequencies 
    % Inputs:
    %   img - ndim filtered image 
    %   ax - filtering axis 
    %   sigma - filtering intensity [0-1 range],  sigma <= 0 no filtering 
    
    
    if sigma <= 0; return ; end
    
    Ndims = ndims(img);
    Npix = size(img);
    shape = ones(1,Ndims);
    shape(ax) = Npix(ax);
    isReal = isreal(img); 
    
    img = fft(img,[],ax);


    x = reshape((-Npix(ax)/2:Npix(ax)/2-1)/Npix(ax), shape);
    spectral_filter = fftshift(exp(1./(-(x.^2)/(sigma^2))));
    img = bsxfun(@times, img, spectral_filter);

    img = ifft(img,[],ax);
    
    if isReal
        img = real(img);
    end


end
