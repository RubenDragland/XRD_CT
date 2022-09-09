function [img, gamma,gamma_x, gamma_y] = stabilize_phase(img, img_ref, mask)
% -----------------------------------------------------------------------
% This file is part of the PTYCHOMAT Toolbox
% Author: Michal Odstrcil, 2016
% License: Open Source under GPLv3
% Contact: ptychomat@gmail.com
% Website: https://bitbucket.org/michalodstrcil/ptychomat
% -----------------------------------------------------------------------
% Description:  adjust phase of the stack of 2D images to be mostly around zero and
% remove linear ramp 

% img - complex image to stabilize
% img_ref - complex images used as reference, if empty ones is used 
% mask  - bool array denoting the region of interest 



    [M,N,~] = size(img);

    if nargin < 2 || isempty(img_ref)
        img_ref = ones(M,N,class(img));
    end
    if nargin < 3
       mask = abs(img_ref) > 0;
    end
     
    
    %% prefer information from the high intensity regions 
    x_ref = bsxfun(@times, img_ref, abs(img_ref));
    x = bsxfun(@times, img , abs(img) );

    %% calculate the optimal phase shift 
    phase_diff = x_ref .* conj(x);
    gamma = mean2(phase_diff .*mask) ./ mean2(mask);
    gamma = gamma./abs(gamma);
    if isnan(gamma); gamma = 1; end
    

    
    %% remove linear phase ramp

    %% initial guess based on  CoM
    [y,x] = center(fftshift_2D(abs(ifft2(phase_diff))).^2, false);
    x=x-floor(M/2)-1;
    y=y-floor(N/2)-1;
    
    xramp = pi*(linspace(-1,1,M));
    yramp = pi*(linspace(-1,1,N));
    
    
    factor = -angle(gamma) + bsxfun(@times,xramp',x)+bsxfun(@times,yramp ,y);
    phase_diff = bsxfun(@times, phase_diff,exp(1i*factor));


    
%%     linear refinement 
    phase_diff= bsxfun(@times, angle(phase_diff), mask); % linearize the problem 
    gamma_x=sum2(bsxfun(@times, (phase_diff),xramp')) ./ sum2(abs(bsxfun(@times,mask,xramp')).^2);
    gamma_y=sum2(bsxfun(@times, (phase_diff),yramp)) ./  sum2(abs(bsxfun(@times,mask,yramp)).^2);
    gamma_x = gamma_x - x;
    gamma_y = gamma_y - y;

    % correct output  image 
    factor = angle(gamma) + bsxfun(@times,xramp',gamma_x)+bsxfun(@times,yramp ,gamma_y);
    img = bsxfun(@times,img,exp(1i*factor)); 

end
