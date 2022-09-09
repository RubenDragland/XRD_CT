function img = imshift_fft(img, x,y, apply_fft)
    % IMSHIFT_FFT  will apply subpixel shift that can be different for 
    %     each frame. If apply_fft == false, then images will be 
    %     assumed to be in fourier space 
    % Author: Michal Odstrcil, 2016
    
    if nargin  < 3
        y = x(:,2);
        x = x(:,1);
    end
    if nargin < 4
        apply_fft = true;
    end
    
    if all(x==0) && all(y==0)
        return
    end
    real_img = isreal(img);
    Np = size(img);

    if apply_fft
         img = fft2(img);
    end

    
    
    xgrid = (fftshift((0:Np(2)-1)/Np(2))-0.5);
    X = reshape((x(:)*xgrid)',1,Np(2),[]);
    X =  exp((-2i*pi)*X);
    img = bsxfun(@times, img,X);
    ygrid = (fftshift((0:Np(1)-1)/Np(1))-0.5);
    Y = reshape((y(:)*ygrid)',Np(1),1,[]);
    Y =  exp((-2i*pi)*Y);
    img = bsxfun(@times, img,Y);
        
    if apply_fft
        img = ifft2(img);
    end
    if real_img
        img = real(img);
    end
    
  
  
end
