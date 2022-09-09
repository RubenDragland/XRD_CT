function [dX, dY] = get_img_grad(img, split, axis)
    %% get vertical and horizontal gradient of the image 
    if nargin < 2 || isempty(split)
        split = 1;
    end
    isReal = isreal(img); 
    Np = size(img);
    if nargin < 3 || any(axis == 1)
        X = 2i*pi*(fftshift((0:Np(2)-1)/Np(2))-0.5);
        dX = bsxfun(@times,fft(img,[],2),X);
        dX = ifft(dX,[],2);
        if isReal; dX = real(dX);end 
    end
    if nargout == 2 || (nargin > 2 && any(axis == 2))
        Y = 2i*pi*(fftshift((0:Np(1)-1)/Np(1))-0.5);
        dY = bsxfun(@times, fft(img,[],1),Y.');
        dY = ifft(dY,[],1);
        if isReal; dY = real(dY);end 
        if nargout == 1; dX = dY; end
    end
end

