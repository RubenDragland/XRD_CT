function img = imrescale_frft(img, scale_x, scale_y)
% -----------------------------------------------------------------------
% This file is part of the PTYCHOMAT Toolbox
% Author: Michal Odstrcil, 2016
% License: Open Source under GPLv3
% Contact: ptychomat@gmail.com
% Website: https://bitbucket.org/michalodstrcil/ptychomat
% -----------------------------------------------------------------------

% IMRESCALE_FRFT based on fractional fourier transformation 
% subpixel precision image rescaling 
    
    isReal = isreal(img);
    if nargin > 2
                
        if any(scale_y ~= 1)
            img = fftshift(ifft(fftshift(FRFT_1D(img,scale_y))));
        end
        if any(scale_x ~= 1)
            img = permute(img,[2,1,3]);
            img = fftshift(ifft(fftshift(FRFT_1D(img,scale_x))));
            img = permute(img,[2,1,3]);
        end
    else
        % 2d version is faster only for many stacked pictures 
        img = fftshift(ifft2(fftshift(FRFT_2D(img,scale_x))));
    end
    if isReal
        img = real(img);
    end
end

function X=FRFT_1D(X,alpha)
    % 1D fractional fourier transformation 
    % See A. Averbuch, "Fast and Accurate Polar Fourier Transform"

    %% it works as magnification lens Claus, D., & Rodenburg, J. M. (2015). Pixel size adjustment in coherent diffractive imaging within the Rayleigh–Sommerfeld regime
    
    %% test plot(abs(fftshift(ifft((FRFT_1D(x,scale))))))
       
    N = size(X,1);
    grid = fftshift(-N:N-1)';

    preFactor = reshape(exp(1i*pi*grid*alpha(:)'),2*N,1,[]);  % perform shift 
    Factor=     reshape(exp(-1i*pi*grid.^2/N * alpha(:)'),2*N,1,[]);  % propagation / scaling
    X=[X; zeros(size(X), class(X))];  % add oversampling 
    X= bsxfun(@times, X, Factor .* preFactor);
    
    % avoid duplication of XX 
    X=fft(X); 
    X = bsxfun(@times, X,fft(conj(Factor)));
    X=ifft(X);
    
    X=bsxfun(@times, X,reshape(Factor .* preFactor,2*N,1,[]));
    X=X(1:N,:,:);
    %% remove phase offset 
    X = bsxfun(@times, X , reshape(exp(-1i*pi*N*alpha/2),1,1,[]));
end

function Y=FRFT_2D(X,alpha)
    % 2D fractional fourier transformation 
    % See A. Averbuch, "Fast and Accurate Polar Fourier Transform"

    %% it maybe works as magification lens Claus, D., & Rodenburg, J. M. (2015). Pixel size adjustment in coherent diffractive imaging within the Rayleigh–Sommerfeld regime

       
    N = size(X,1);
    grid = fftshift(-N:N-1)';
    [Xg,Yg] = meshgrid(grid, grid);
    preFactor = exp(1i*pi*(Xg+Yg)*alpha);  % perform shift after FFT 
    Factor=exp(-1i*pi*(Xg.^2+Yg.^2)/N(1) * alpha);  % propagation / scaling
    x_tilde = zeros(2*N, 2*N, size(X,3),  class(X));
    x_tilde(1:N, 1:N,:) = X;
    x_tilde= bsxfun(@times, x_tilde, Factor .* preFactor);
        
    XX=fft2(x_tilde);

    YY=fft2(conj(Factor));  

    XX = bsxfun(@times, XX,YY);

    Y=ifft2( XX );
    Y=bsxfun(@times, Y,Factor .* preFactor );
    Y=Y(1:N,1:N,:);

end

