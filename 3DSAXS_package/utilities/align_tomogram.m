function  delta =  align_tomogram(tomogram, projection, which_one);
%*-------------------------------------------------------------------------------------*
%|                                                                                                           |
%|  Except where otherwise noted, this work is licensed under a            |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0        |
%|  International (CC BY-NC-SA 4.0) license.                                         |
%|                                                                                                           |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)       |
%|                                                                                                           |
%|      Author: CXS group, PSI                                                                |
%*------------------------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form and additionally you should cite:
%     M. Liebi, M. Georgiadis, A. Menzel, P. Schneider, J. Kohlbrecher, 
%     O. Bunk, and M. Guizar-Sicairos, “Nanostructure surveys of 
%     macroscopic specimens by small-angle scattering tensor tomography,”
%     Nature 527, 349-352 (2015).   (doi:10.1038/nature16056)
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
   % Size of tomogram
    N = size(tomogram);
    % defines the center in each direction
    % X appears in the second position of the meshgrid,
    % because meshgrid is (x,y,z) and size is (y,x,z)
    x = [1:N(2)] - ceil(N(2)/2); %
    y = [1:N(1)] - ceil(N(1)/2);
    z = [1:N(3)] - ceil(N(3)/2);
    
    % creates a mesh for arb_projection function
    [X, Y, Z] = meshgrid(x,y,z);
    
  
    tosee = mean(projection(which_one).data, 3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p.volume_upsampling = 1;  % Upsamples the volume by a factor of 2 to reduce artifacts
    p.method = 'bilinear';    % 'nearest' or 'bilinear' method for volume allocation
    p.filter_2D = 3;          % Strength of filter applied to image [0-3]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Register the measured projection and the projection from the tomogram
    xout = [1:size(tosee,2)]- ceil(size(tosee,2)/2);
    yout = [1:size(tosee,1)]- ceil(size(tosee,1)/2);
    % calculate the reference image based on the rotation matrix Rot_exp_now
    [proj_out_all, xout, yout] = arb_projection(tomogram,X,Y,Z,...
        projection(which_one).Rot_exp, p,xout,yout);
    
    % reference image generated from the tomogram in STEP 2
    ref_image =  proj_out_all(:,:);
    % normalize by the maximum
    %ref_image = ref_image./(max(max(ref_image)));

    % define the images based on the output of diode or sym_int data
    to_be_aligned = tosee;
   
   
    upsamp = 100;
    % displ     = 0  no information displayed (default)
    %           = 1  to display text information
    %           > 1  also shows images in figure(display)
    displ = 0;
    % Wfilt:     Fourier domain filter for registration (= 1 for no effect)
    Wfilt = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fromx = [];
tox = [];
fromy = [];
toy = [];
    [subim1, subim2, delta, ~,~] = registersubimages_2...
        (ref_image, to_be_aligned, 1:length(xout), fromy:toy, ...
        fromx:tox, fromy:toy, upsamp, displ, Wfilt);