function img = smooth_edges(img, smooth_kernel)
    % Function: img = smooth_edges(img, smooth_kernel)
    %   takes stack of 2D images and smooths boundaries to avoid sharp edge artefacts during imshift_fft 
    %   Inputs:
    %       img - 2D / 3D array, smoothing is done along first two dimensions 
    %       smooth_kernel - size of the smoothing region, default is 3 
    
    if nargin < 2
        smooth_kernel = 3; 
    end
    
    Npix = size(img);
    for i = 1:2
        img =  circshift(img,ceil(Npix(i)/2),i); 
        ind = {':',':',':'};
        ind{i} = ceil(Npix(i)/2)+(-smooth_kernel:smooth_kernel); 
        ker_size = [1,1]; 
        ker_size(i) = smooth_kernel;
        img_tmp = img(ind{:}); 
        %% smooth across the image edges 
        img_tmp = convn(img_tmp, ones(ker_size) , 'same');
        %% avoid boundary issues 
        boundary_shape = [1,1]; 
        boundary_shape(i) = 2*smooth_kernel+1; 
        img_tmp = bsxfun(@rdivide, img_tmp, conv(ones(boundary_shape), ones(ker_size), 'same'));
        img(ind{:}) = img_tmp; 
        img =  circshift(img,ceil(Npix(i)/2),i); 
    end

end