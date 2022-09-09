function img = imshift_fft_ax(img, shift, ax, apply_fft)
% -----------------------------------------------------------------------
% This file is part of the PTYCHOMAT Toolbox
% Author: Michal Odstrcil, 2016
% License: Open Source under GPLv3
% Contact: ptychomat@gmail.com
% Website: https://bitbucket.org/michalodstrcil/ptychomat
% -----------------------------------------------------------------------

% IMSHIFT_FFT_AX  will apply subpixel shift that can be different
%     for each frame along axis AX. If apply_fft == false, 
%     then images will be assumed to be in fourier space 



    if nargin < 4
        apply_fft = true;
    end
    if all(shift == 0)
        return
    end
    
    isReal = isreal(img); 
    
    Npix = size(img);
    
    if ndims(img) == 3
        Np = [1,1,Npix(3)];
    else
        Np = Npix;
        Np(ax) = 1; 
    end

    Ng = ones(1,3);
    Ng(ax) = Npix(ax);
    
    
    
    if isscalar(shift)
         shift = shift .* ones(Np);
    end
    

    grid = single(fftshift((0:Npix(ax)-1)/Npix(ax))-0.5);
    
    if isa(img, 'gpuArray')
        grid = gpuArray(grid);
    end
    
    X = bsxfun(@times, reshape(shift,Np), reshape(grid,Ng));
    X =  exp((-2i*pi)*X);
    
    
    if apply_fft
        img = fft(img,[],ax);
    end
        
    img = bsxfun(@times, img,X);
   
    if apply_fft
        img = ifft(img,[],ax);
    end
        
    if isReal
        img = real(img);
    end

end

