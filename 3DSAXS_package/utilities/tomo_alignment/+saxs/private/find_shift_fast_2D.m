function shift = find_shift_fast_2D(o1, o2, sigma, apply_fft)   
    
    % FIND_SHIFT_FAST_2D uses crosscorelation to find shift between o1 nd
    % o2 patterns 
    % Inputs:
    %   o1 - aligned array 2D or 3D
    %   o2 - template for alignement 2D or 3D
    %   sigma - filtering intensity [0-1 range],  sigma <= 0 no filtering 
    %   apply_fft - if false, assume o1 and o2 to be already in fourier domain 
    
    
   if nargin < 4 
       apply_fft = true;
   end
   
   if apply_fft

        if size(o2,3) == 1
            % avoid issues caused by empty space around 
            [~,ROI] = get_ROI(o2>0,-0.05);
            o1 = o1(ROI{:},:);
            o2 = o2(ROI{:},:);
        end


        [nx, ny, ~] = size(o1);

        % suppress edge effects of the registration procedure
        spatial_filter = tukeywin(nx) * tukeywin(ny)'; 

        o1 = bsxfun(@times, o1, spatial_filter);
        o2 = bsxfun(@times, o2, spatial_filter);

        o1 = fft2(o1);
        o2 = fft2(o2);
   end
   
    [nx, ny, ~] = size(o1);

    if sigma > 0 
        % remove low frequencies 
        [X,Y] = meshgrid( (-nx/2:nx/2-1), (-ny/2:ny/2-1));
        spectral_filter = fftshift(exp(1./(-(X.^2+Y.^2)/(mean([nx,ny])*sigma)^2)))';
        o1 = bsxfun(@times, o1, spectral_filter);
        o2 = bsxfun(@times, o2, spectral_filter);       
    end
    
    % fast subpixel cross correlation 
    xcorrmat = fftshift(abs(ifft2(o1.*conj(o2)))); 
    mxcorr = mean(xcorrmat,3); 
    [m,n] = find(mxcorr == max(mxcorr(:)));

    

    MAX_SHIFT = 10;  % +-10px search 
    idx = {[m-MAX_SHIFT:m+MAX_SHIFT],[n-MAX_SHIFT:n+MAX_SHIFT],':'};
    xcorrmat = xcorrmat(idx{:});
    MAX = max(max(xcorrmat));
    
    xcorrmat = bsxfun(@times, xcorrmat, 1. / MAX);

    %% get CoM of  the central  peak only !!, assume a single peak 
    xcorrmat = max(0, xcorrmat - 0.5).^4;
    [x,y] = find_center_fast(xcorrmat);
    shift = [x,y]+[n,m]-floor([ny,nx]/2)-1;
end

function [x,y,MASS] = find_center_fast(xcorrmat)
    MASS = squeeze(sum(sum(xcorrmat)));
    [N,M,~] = size(xcorrmat); 
    x = squeeze(sum( bsxfun(@times, sum(xcorrmat,1), 1:M), 2)) ./ MASS - floor(M/2)-1;
    y = squeeze(sum(bsxfun(@times, sum(xcorrmat,2), (1:N)'),1)) ./ MASS - floor(N/2)-1;
end
