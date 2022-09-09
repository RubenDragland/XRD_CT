function [ proj_out_all, xout, yout ] = arb_projection( tomo_obj_all, X, Y, Z, R, p, xout, yout )
% [ proj_out_all, xout, yout ] = arb_projection( tomo_obj_all, X, Y, Z, R, p, xout, yout )
%
% Computes a projection along Z after an arbitrary rotation of the object
% tomo_obj_all defined by the rotation matrix R.
%
% Inputs
%  tomo_obj_all   - A 3D or 4D array to allow for tensor tomography, first
%                   three indices are spatial coordinates and the fourth allows for instance
%                   for another quantity that varies along the fourth index, e.g. detector
%                   angular section.
%  X, Y, Z        - 3D arrays that define the array dimension that
%                   corresponds to X Y Z respectively. They must be integer
%                   valued and they must have unit interval between voxels.
%                   May be created with meshgrid. These arrays allow you to
%                   define which array dimension corresponds to X, Y, Z
%                   respectively. The rotation is around X = Y = Z = 0.
%                   Even if this point does not exist in the array.
% R               - 3x3 rotation matrix that describes the rotation of the
%                   object in vector convention. [x;y;z] = R*[x;y;z]
% p               - Structure that defines the parameters of the projection
% p.volume_upsampling = 1 - the volume is upsampled by a factor 2 before
%                   the projection. Algorithm takes about 8 times longer
%                   but reduces artifacts significantly without sacrificing
%                   too much accuracy. Default = 0.
% p.filter_2D       Filter strength on the resulting 2D projection to reduce
%                   artifacts but also resolution. ( = 0) No filter.
%                   Maximum filtering for ( = 3) default. Note that if you
%                   use p.volume_upsampling = 1 the filter is applied on an
%                   upsampled 2D projection, so the decrease of resolution
%                   by the filter will be much less.
% p.method          'nearest' or 'bilinear', the method affects the
%                   allocation of the voxel values into the 2D array.
% xout, yout        Row arrays, e.g. [-10:10] that specify the output
%                   coordinates of the projection. Values should be spaced
%                   by one but a subpixel shift can be added, could be useful
%                   for including projection misalignment. By default
%                   the vectors corresponds to the range of X and Y.
%
% Manuel Guizar - 2014.06.05
%
% If not using recommended settings significant artifacts have been observed
% for rotations by 23, 28, 37, 45 or 53 degrees. I recommend to look once
% at a projection of a sphere from the angles of your experiment to
% determine how strong the artifacts would be. Code to generate a smooth
% sphere is included as comments at the end of the function.
%
% Recommended settings:
%%%% Fast & low artifacts & low resolution %%%
%   p.volume_upsampling = 0;
%   p.filter_2D = 3;           % Lower values give you compromise between resolution and artifacts
%   p.method = 'bilinear';     % 'nearest' or 'bilinear'
%
%%%% Slow & low artifacts & high resolution %%%
%   Still gives some artifacts at 45 degrees if only rotated around z
%   p.volume_upsampling = 1;
%   p.filter_2D = 3;           % Lower values give you compromise between resolution and artifacts
%   p.method = 'bilinear';     % 'nearest' or 'bilinear'

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

if isfield(p,'volume_upsampling')
    volume_upsampling = p.volume_upsampling;  % Upsamples the volume by a factor of 2 to reduce artifacts
else
    volume_upsampling = 0;
end

if isfield(p,'method')
    method = p.method;
else
    method = 'bilinear';     % 'nearest' or 'bilinear' method for volume allocation
end

if isfield(p,'filter_2D')
    filter_2D = p.filter_2D;
else
    filter_2D = 3; % Strength of filter applied to image
end

% Optional input window, there can be a subpixel offset
% if nargin == 6
%     xout = [-20:20];
%     yout = [-30:30];
% end

if (~all(size(tomo_obj_all(:,:,:,1))==size(X)))||(~all(size(tomo_obj_all(:,:,:,1))==size(Y)))||(~all(size(tomo_obj_all(:,:,:,1))==size(Z)))
    error('Size of Tomogram and input X, Y, Z does not match')
end

% Check that spacing between voxels is 1 and that they are oriented along
% main matrix axes. This is an assumption taken for the bilinear
% interpolation
if ~((abs(X(2,2,2)-X(1,1,1))==1)&&(abs(Y(2,2,2)-Y(1,1,1))==1)&&(abs(Z(2,2,2)-Z(1,1,1))==1))
    error('Spacing between consecutive coordinates is not 1')
end

if ~exist('xout','var')||~exist('yout','var')
    % Determine directions along which X and Y run
    if abs(X(2,1,1)-X(1,1,1))==1
        xdir = 1;
    elseif abs(X(1,2,1)-X(1,1,1))==1
        xdir = 2;
    elseif abs(X(1,1,2)-X(1,1,1))==1
        xdir = 3;
    else
        error('Could not determine along which array dimension X runs')
    end
    if abs(Y(2,1,1)-Y(1,1,1))==1
        ydir = 1;
    elseif abs(Y(1,2,1)-Y(1,1,1))==1
        ydir = 2;
    elseif abs(Y(1,1,2)-Y(1,1,1))==1
        ydir = 3;
    else
        error('Could not determine along which array dimension Y runs')
    end
    
    % Default definition of output coordinates
    switch xdir
        case 1
            xout = squeeze(X(:,1,1));
        case 2
            xout = squeeze(X(1,:,1));
        case 3
            xout = squeeze(X(1,1,:));
    end
    switch ydir
        case 1
            yout = squeeze(Y(:,1,1));
        case 2
            yout = squeeze(Y(1,:,1));
        case 3
            yout = squeeze(Y(1,1,:));
    end
else
    % if ~((xout(2)-xout(1)==1)&&(yout(2)-yout(1)==1)) %%rounding issues
    if ~(abs(yout(2)-yout(1))-1<1e-5)&&(abs(xout(2)-xout(1))-1<1e-5)
        error('xout and yout must be monotonically increasing in integer steps')
    end
end

% Generate coordinates for rotated object
Xp = R(1,1)*X + R(1,2)*Y + R(1,3)*Z;
Yp = R(2,1)*X + R(2,2)*Y + R(2,3)*Z;
% Zp = R(3,1)*X + R(3,2)*Y + R(3,3)*Z; % Not needed


%%% Here I assign for each 4th index a tomogram and operate element by
%%% element on the 3D array. I think this is the most efficient because the
%%% first three indices are continuous in memory. However it remains to be
%%% explored if there is a way of looping through the 3 indices of the
%%% tomogram and operating on all elements of the fourth index
%%% simultaneously
if volume_upsampling
    % Upsample by factor 2
    Xp = trilinupsamp(2*Xp);
    Yp = trilinupsamp(2*Yp);
    
    xout_orig = xout;
    xout = 2*xout_orig(1) : 2*xout_orig(end);
    yout_orig = yout;
    yout = 2*yout_orig(1) : 2*yout_orig(end);
    
    tomo_obj_all = trilinupsamp(tomo_obj_all);
end

% shortcuts
min_xout = min(xout);
min_yout = min(yout);

if strcmpi(method,'nearest')
    %%% Nearest neighbor %%%
    % Assigns value of voxel to the nearest pixel
    Ax = round(Xp-min_xout+1); % Index of output image that corresponds to the voxel
    Ay = round(Yp-min_yout+1); % Index of output image that corresponds to the voxel
    
    % initialization for all iterations of the following loop
    proj_out_all = zeros(numel(yout), numel(xout), size(tomo_obj_all, 4));

    for jj = 1:size(tomo_obj_all,4)
        for ii = 1:numel(Ax)
            if (Ax(ii) > 0)&&(Ax(ii) < size(proj_out,2))&&(Ay(ii) > 0)&&(Ay(ii) < size(proj_out,1))
                proj_out(Ay(ii),Ax(ii)) = proj_out(Ay(ii),Ax(ii)) + tomo_obj(ii);
            end
        end
        proj_out_all(:,:,jj) = proj_out;
    end
    
elseif strcmpi(method,'bilinear')
    %%% Bilinear interpolation %%%
    Ax = floor(Xp-min_xout+1); % Index of output image that corresponds to the voxel
    Ay = floor(Yp-min_yout+1); % Index of output image that corresponds to the voxel
    Tx = (Xp-min_xout+1) - Ax; % Variable from 0 to 1 from x distance of pixel Ax to Ax+1 where the voxel hits
    Ty = (Yp-min_yout+1) - Ay;
    
    size_proj_out_all = uint64([numel(yout), numel(xout), size(tomo_obj_all, 4)]);
    proj_out_all = sum_projection(size_proj_out_all, tomo_obj_all, Ax, Ay, Tx, Ty);
    
else
    error('Method was not recognized')
end

%%% Adding 2D filtering at the end. The triangular filter emulates
%%% bilinear upsampling and down sampling and seems very effective
switch filter_2D
    case 0
    case 1
        filt_2D = [0.25 1 0.25];
    case 2
        filt_2D = [0.5 1 0.5];
    case 3
        filt_2D = [1/3 2/3 1 2/3 1/3];
    otherwise
        error('Filter 2D specified is not recognized')
end
if filter_2D > 0
    filt_2D = filt_2D.'*filt_2D;
    filt_2D = filt_2D/sum(filt_2D(:));
    for jj = 1:size(tomo_obj_all,4)
        proj_out_all(:,:,jj) = conv2(proj_out_all(:,:,jj), filt_2D, 'same');
    end
end

%%% Downsampling %%%
% Currently done as a nearest neighbor, seems quite good but may be better
% if a bilinear approach was used.
if volume_upsampling
    proj_out_all = proj_out_all(1:2:end,1:2:end,:)/2;
    xout = xout_orig;
    yout = yout_orig;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% How to make a test sphere with smooth edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Define coordinates, center of rotation at zero
% N = [100 50 50]; % Size of tomogram
% x = [1:N(2)]-ceil(N(2)/2);
% y = [1:N(1)]-ceil(N(1)/2);
% z = [1:N(3)]-ceil(N(3)/2);
% [X Y Z] = meshgrid(x,y,z);
%
% % Define object
% % tomo_obj = double((X.^2+Y.^2+Z.^2)<=20^2);
% % tomo_obj = tomo_obj+double((X.^2+(Y-20).^2+Z.^2)<=10^2);
% tomo_obj_all(:,:,:,1) = exp(-((X.^2+Y.^2+Z.^2)./20^2).^20);
% tomo_obj_all(:,:,:,1) = tomo_obj_all(:,:,:,1)+exp(-((X.^2+(Y-20).^2+Z.^2)./10^2).^20);
% tomo_obj_all(:,:,:,2) = exp(-((X.^2+Y.^2+Z.^2)./10^2).^20);
% tomo_obj_all(:,:,:,2) = tomo_obj_all(:,:,:,2)+exp(-(((X-5).^2+(Y+10).^2+Z.^2)./10^2).^20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Other kernels that were tried for the allocation, they did not work
%%% particularly well, the sinc was quite good but slower than the volume
%%% interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % An attempt to generalize the kernel, only written here for a 3x3 kernel
% kern = [1/4 1/2 1/4];
% kern = kern'*kern;
% %%% Kernel interpolation %%%
% Ax = round(Xp-min(xout)+1); % Index of output image that corresponds to the voxel
% Ay = round(Yp-min(yout)+1); % Index of output image that corresponds to the voxel
% % Intialize output array
% proj_out = zeros(numel(yout),numel(xout));
% for ii = 1:numel(Ax)
%     if (Ax(ii) > 1)&&(Ax(ii) < size(proj_out,2))&&(Ay(ii) > 1)&&(Ay(ii) < size(proj_out,1))
%         for jj = 1:size(kern,1)
%             for kk = 1:size(kern,2)
%                 proj_out(Ay(ii)+jj-2,Ax(ii)+kk-2) = proj_out(Ay(ii)+jj-2,Ax(ii)+kk-2) + tomo_obj(ii)*kern(jj,kk);
%             end
%         end
%     end
% end

% % An attempt to generalize the kernel, only written here for a 3x3 kernel
% % Including shift with parabolic interpolator - did not work very well this
% % interpolator
% %%% Kernel interpolation %%%
% Ax = round(Xp-min(xout)+1); % Index of output image that corresponds to the voxel
% Ay = round(Yp-min(yout)+1); % Index of output image that corresponds to the voxel
% Tx = (Xp-min(xout)+1) - Ax;               % Variable from -0.5 to 0.5 from x distance of pixel Ax to Ax+1 where the voxel hits
% Ty = (Yp-min(yout)+1) - Ay;
% % Intialize output array
% proj_out = zeros(numel(yout),numel(xout));
% kernx = [-1 0 1;-1 0 1;-1 0 1]*0.3536;
% kerny = [-1 -1 -1; 0 0 0; 1 1 1]*0.3536;
% for ii = 1:numel(Ax)
%     if (Ax(ii) > 1)&&(Ax(ii) < size(proj_out,2))&&(Ay(ii) > 1)&&(Ay(ii) < size(proj_out,1))
%         % This is the sinc^2 kernel, maybe I could try with a parabolic
%         % approximation of this? Could be faster and yield less artifacts
%         % Still shows some artifacts for special angles
%         kern = 2-2*((kernx-Tx(ii)).^2+(kerny-Ty(ii)).^2).^3;
%         kern = kern/sum(kern(:));
%         for jj = 1:size(kern,1)
%             for kk = 1:size(kern,2)
%                 proj_out(Ay(ii)+jj-2,Ax(ii)+kk-2) = proj_out(Ay(ii)+jj-2,Ax(ii)+kk-2) + tomo_obj(ii)*kern(jj,kk);
%             end
%         end
%     end
% end

% % An attempt to generalize the kernel, only written here for a 3x3 kernel
% % Including shift with sinc^2 interpolator
% % Works well, but its so slow
% %%% Kernel interpolation %%%
% Ax = round(Xp-min(xout)+1); % Index of output image that corresponds to the voxel
% Ay = round(Yp-min(yout)+1); % Index of output image that corresponds to the voxel
% Tx = (Xp-min(xout)+1) - Ax;               % Variable from -0.5 to 0.5 from x distance of pixel Ax to Ax+1 where the voxel hits
% Ty = (Yp-min(yout)+1) - Ay;
% % Intialize output array
% proj_out = zeros(numel(yout),numel(xout));
% for ii = 1:numel(Ax)
%     if (Ax(ii) > 1)&&(Ax(ii) < size(proj_out,2))&&(Ay(ii) > 1)&&(Ay(ii) < size(proj_out,1))
%         % This is the sinc^2 kernel, maybe I could try with a parabolic
%         % approximation of this? Could be faster and yield less artifacts
%         % Still shows some artifacts for special angles
%         kernx = [-1/2 0 1/2]-Tx(ii)/2+rand*1e-6;
%         kernx = (sin(pi*kernx)./(pi*kernx)).^2;
%         kerny = [-1/2 0 1/2]-Ty(ii)/2+rand*1e-6;
%         kerny = (sin(pi*kerny)./(pi*kerny)).^2;
%         kern = kerny'*kernx;
%         kern = kern/sum(kern(:));
%         for jj = 1:size(kern,1)
%             for kk = 1:size(kern,2)
%                 proj_out(Ay(ii)+jj-2,Ax(ii)+kk-2) = proj_out(Ay(ii)+jj-2,Ax(ii)+kk-2) + tomo_obj(ii)*kern(jj,kk);
%             end
%         end
%     end
% end
% % Best results jet but it takes 5 seconds, 14 times longer than nearest
% % neighbor and still gives artifacts



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Good performance was obtained for downsampling followed by upsampling
%%% on the resulting projection.
%%% Below I explored that this is equivalent to a triangular filter on the
%%% output image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Playing around with downsample, then upsample
% downsamplemethod = 'bilinear';
% % Downsample
% downsampling = 3;
% % proj_out_orig = proj_out;
% % proj_out = proj_out_orig(1:downsampling:end-downsampling+1,1:downsampling:end-downsampling+1);
% % for ii = 1:downsampling
% proj_out = imresize(proj_out,1/downsampling,downsamplemethod);
% proj_out = imresize(proj_out,downsampling,downsamplemethod);
%
% auxfilt = zeros(20);
% auxfilt(11,11) = 1;
% auxfilt = imresize(auxfilt,1/downsampling,downsamplemethod);
% auxfilt = imresize(auxfilt,downsampling,downsamplemethod);
% figure; imagesc(auxfilt); axis xy equal tight


% function [proj_out] = sum_projection(proj_out,tmp_3D,A_linear_index)
%
% for ind=1:size(A_linear_index,1)
%     proj_out(A_linear_index(ind)) = proj_out(A_linear_index(ind)) + tmp_3D(ind);
% end
%
%
