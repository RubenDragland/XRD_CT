%HSV_VOLUMES_2RGB Function for creating colorful pictures of volumes.
%    Useful for plotting colocalized values. It receives three volumes that
%    represent the local Hue, Saturation and Value. Based on the given
%    parameters, chosen color range and thresholds a volume of RGB values
%    is returned.
%
% [RGBvol, imhan, cbarhan] = hsv_volumes_2rgb(H,S,V,params)
%
% Inputs:
% H, S, V       Volumes to represent the hue, saturation, and value
%               respectively
%
% *optional*
% ** params     Structure with parameters
% **   params.value_thresh  Min and max values assigned to Value (or Intensity) of 0 and 1 respectively. eg. = [0.2 0.9]. Default are the minimun and maximum values in the V volume.
% **   params.sat_thresh    Min and max values assigned to Saturation of 0 and 1 respectively. eg. = [0 1.3]. Default are the minimun and maximum values in the S volume.
% **   params.hue_thresh    Min and max values assigned to minimum and maximum Hue values, respectively. eg. = [5.477 5.48]. Default are the minimun and maximum values in the H volume. More detailed explanation is given below.
% **   params.hue_val       Min and max hue used in the plot. Default = [0 1];
% **   params.plotting_flag Plots the result. Default = false.
% **   params.init_frame    Choose the frame shown in the plot. Default is the middle of the stack.
% **   params.image_filename  Path and filename string to save images of the figures shown. Default = []. Eg. lemon_%s.png. Give a %s in the name since there are two figures generated and this is used to distinguish them.
% **   params.bmp_filename Path and filename to save a color BMP stack.Default = []. e.g '/das/work/units/csaxs/p17880/TIFF/pristine_lemon/lemon_%05d.bmp'
%
% returns
% ++ RGBvol     Volume with RGB values
% ++ imhan      Figure handle for the image figure
% ++ cbarhan    Figure handle for the 3D colorbar figure
%
%
% see also: plotting.imagesc3D
%
% EXAMPLES:
%   params.value_thresh = [0.2 0.9];         % Which values will have a intensity (value) of zero and one
%   params.hue_thresh   = [4.5 6.2];  % It will represent values of H from 4.5 to 6.2
%   params.hue_val      = [0 0.5];   % This controls the hues that will be displayed, in this case 0 (red) to 0.5 (cyan).
%   params.sat_thresh    = [0 1.3];
%   params.plotting_flag    = true;
%   params.image_filename   = '/path/images/orange%s.png';
%   params.tiff_filename    = '/path/pristine_orange/orange.tiff';
%
%   [RGBvol, imhan, cbarhan] = hsv_volumes_2rgb(H,S,V,params)


%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2018 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%*-----------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing language this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS software package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.

function [RGBvol, imhan, cbarhan] = hsv_volumes_2rgb(H,S,V,params)

value_vol = V;
hue_vol   = H;
sat_vol   = S;

if isfield(params,'value_thres')
    value_thresh = params.value_thresh;
else
    value_thresh = [min(V(:)) max(V(:))];
end

if isfield(params,'sat_thresh')
    sat_thresh = params.sat_thresh;
else
    sat_thresh = [min(S(:)) max(S(:))];
end

if isfield(params,'hue_thresh')
    hue_thresh = params.hue_thresh;
else
    hue_thresh = [min(H(:)) max(H(:))];
end

if isfield(params,'hue_val')
    hue_val = params.hue_val;
else
    hue_val = [0 1];
end

if isfield(params,'plotting_flag')
    plotting_flag = params.plotting_flag;
else
    plotting_flag = false;
end

if isfield(params,'image_filename')
    image_filename = params.image_filename;
    if ~plotting_flag
        error('Since you specify params.png_filename you should set the params.plotting_flag to true')
    end
else
    image_filename = [];
end

if isfield(params,'bmp_filename')
    bmp_filename = params.bmp_filename;
else
    bmp_filename = [];
end

if isfield(params,'init_frame')
    init_frame = params.init_frame;
else
    init_frame = round(size(value_vol,3)/2);
end

% Apply threshold and normalize
value_vol_thresh = value_vol;   
value_vol_thresh(value_vol_thresh<value_thresh(1)) = value_thresh(1);
value_vol_thresh(value_vol_thresh>value_thresh(2)) = value_thresh(2);
value_vol_thresh = value_vol_thresh - min(value_vol_thresh(:));
value_vol_thresh = value_vol_thresh./max(value_vol_thresh(:));

hue_vol_thresh = hue_vol;
hue_vol_thresh(hue_vol_thresh<hue_thresh(1)) = hue_thresh(1);
hue_vol_thresh(hue_vol_thresh>hue_thresh(2)) = hue_thresh(2);
hue_vol_thresh = hue_vol_thresh - min(hue_vol_thresh(:));
hue_vol_thresh = hue_vol_thresh./max(hue_vol_thresh(:));
hue_vol_thresh = hue_vol_thresh.*diff(hue_val) + hue_val(1);

sat_vol_thresh = sat_vol;
sat_vol_thresh(sat_vol_thresh<sat_thresh(1)) = sat_thresh(1);
sat_vol_thresh(sat_vol_thresh>sat_thresh(2)) = sat_thresh(2);
sat_vol_thresh = sat_vol_thresh - min(sat_vol_thresh(:));
sat_vol_thresh = sat_vol_thresh./max(sat_vol_thresh(:));

% Colorbar


hsv_vol = cat(4,hue_vol_thresh,sat_vol_thresh,value_vol_thresh);
hsv_vol(:,:,:,1) = mod(hsv_vol(:,:,:,1),1);
N = size(hsv_vol);
hsv_vol = reshape(hsv_vol,[N(1)*N(2)*N(3) 3]);
RGBvol = hsv2rgb(hsv_vol);
RGBvol = reshape(RGBvol,N);


% Plotting
if plotting_flag
    
    hsv_vol_color = hsv_vol;
    hsv_vol_color(:,2:3) = 1;
    rgb_vol_color = hsv2rgb(hsv_vol_color);
    rgb_vol_color = reshape(rgb_vol_color,N);

    imhan = figure(4);
    clf
    ax(1) = subplot(2,2,1);
    plotting.imagesc3D(value_vol,'init_frame',init_frame)
    colormap(bone)
    axis xy equal tight
    colorbar
    caxis(value_thresh)
    title('Normalized Thresholded Value')
    
    ax(2) = subplot(2,2,2);
    plotting.imagesc3D(rgb_vol_color,'init_frame',init_frame)
    axis xy equal tight
    title('Normalized Thresholded Hue')
    
    ax(3) = subplot(2,2,3);
    % plotting.imagesc3D(sat_vol_thresh,'init_frame',init_frame)
    plotting.imagesc3D(sat_vol,'init_frame',init_frame); caxis(sat_thresh)
    colormap(bone)
    axis xy equal tight
    colorbar
    title('Normalized Thresholded Saturation')
    
    ax(4) = subplot(2,2,4);
    plotting.imagesc3D(RGBvol,'init_frame',init_frame)
    axis xy equal tight
    title('Colorful')
    
    cmap_ax = axes('Position',[0.93 0.6 0.015 0.31],'XTick',[]);
    mapa = zeros(64,2,3);
    mapa(:,1,1) = linspace(hue_val(1),hue_val(2),64);
    mapa(:,1,1) = mod(mapa(:,1),1);
    mapa(:,1,2) = 1;
    mapa(:,1,3) = 1;
    mapa(:,2,:) = mapa(:,1,:);
    mapa = hsv2rgb(mapa);
    imagesc([1 2],hue_thresh,mapa);
    set(cmap_ax,'XTick',[]);
    axis xy
    
    % Plotting 3D colorbar
    cbarhan = figure(2);
    clf
    axh = subplot(1,3,1);%axes;
    % First, create a 100-by-100 image to texture the cone with:
    hold on
    % H = repmat(linspace(0, 1, 100), 100, 1);     % 100-by-100 hues
    H = repmat(linspace(hue_val(1), hue_val(2), 100), 100, 1);     % 100-by-100 hues
    S = repmat([linspace(0, 1, 50) ...           % 100-by-100 saturations
        linspace(1, 0, 50)].', 1, 100);  %'
    V = repmat([ones(1, 50) ...                  % 100-by-100 values
        linspace(1, 0, 50)].', 1, 100);  %'
    hsvImage = cat(3, H, S, V);                  % Create an HSV image
    C = hsv2rgb(hsvImage);                       % Convert it to an RGB image
    
    % Next, create the conical surface coordinates:
    % theta = linspace(0, 2*pi, 100);  % Angular points
    theta = linspace(0, 3*pi/2, 100);  % Angular points
    X = [zeros(1, 100); ...          % X coordinates
        cos(theta); ...
        zeros(1, 100)];
    Y = [zeros(1, 100); ...          % Y coordinates
        sin(theta); ...
        zeros(1, 100)];
    Z = [2.*ones(2, 100); ...        % Z coordinates
        zeros(1, 100)];
    
    % Finally, plot the texture-mapped surface:
    surf(X, Y, Z, C, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
    % surf(X(1:3,:), Y(1:3,:), Z(1:3,:),  'FaceColor', 'texturemap', 'EdgeColor', 'none');
    axis equal
    
    
    % Create the coordinates for the flat sides of the cut cone
    Xt = [zeros(1, 100); ...          % X coordinates
        linspace(0,1,100); ...
        zeros(1, 100)];
    Yt = [zeros(1, 100); ...          % Y coordinates
        zeros(1, 100); ...
        zeros(1, 100)];
    Zt = [2.*ones(2, 100); ...        % Z coordinates
        zeros(1, 100)];
    
    St = repmat(linspace(0, 1, 100), 100, 1);
    
    hsvImaget = cat(3, H*0+hue_val(1), St, V);                  % Create an HSV image
    Ct = hsv2rgb(hsvImaget);                       % Convert it to an RGB image
    hsvImagett = cat(3, H*0+hue_val(2), St, V);                  % Create an HSV image
    Ctt = hsv2rgb(hsvImagett);                       % Convert it to an RGB image
    
    surf(Xt, Yt, Zt,Ct,  'FaceColor', 'texturemap', 'EdgeColor', 'none');
    
    surf(Yt, -Xt, Zt,Ctt,  'FaceColor', 'texturemap', 'EdgeColor', 'none');
    axis equal
    axh.View = [13.6000 18.8000];
    hold off
    
    axh2 = subplot(1,3,2);
    surf(Xt, Yt, Zt,Ct,  'FaceColor', 'texturemap', 'EdgeColor', 'none');
    axis equal
    axh2.View = [0 0];
    
    axh3 = subplot(1,3,3);
    surf(Yt, -Xt, Zt,Ctt,  'FaceColor', 'texturemap', 'EdgeColor', 'none');
    axis equal
    axh3.View = [-90 0];
end

if ~isempty(image_filename)
    figure(imhan);
    savename = sprintf(image_filename,'_image.png');
    fprintf('Saving %s\n',savename);
    saveas(gcf,savename);
    figure(cbarhan);
    savename = sprintf(image_filename,'_colorbar.png');
    saveas(gcf,savename);
    fprintf('Saving %s\n',savename);
end

if ~isempty(bmp_filename)
    
    [filepath, ~, ~] = fileparts(bmp_filename);
    
    if ~exist(filepath,'dir')
        mkdir(filepath)
    end
       
    for ii = 1:size(RGBvol,3)
        fprintf('Writting %s\n', sprintf(bmp_filename,ii))
        imwrite(squeeze(RGBvol(:,:,ii,:)),sprintf(bmp_filename,ii))
    end

    %     including ImageWidth, ImageLength, BitsPerSample, 
%     SamplesPerPixel, Compression, PlanarConfiguration, and Photometric.

%     obj = Tiff(tiff_filename,'w8');
%     for ii = 1:size(RGBvol,3)
% %         obj.nextDirectory();
%         obj.writeDirectory();
%         setTag(obj,'ImageLength',size(RGBvol,1));
%         setTag(obj,'ImageWidth' ,size(RGBvol,2));
%         setTag(obj,'Photometric',Tiff.Photometric.RGB);
%         setTag(obj,'BitsPerSample',16);
%         setTag(obj,'SamplesPerPixel',3);
%         setTag(obj,'Compression',Tiff.Compression.None);
%         setTag(obj,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
%         obj.write(squeeze(uint16(RGBvol(:,:,ii,:))));
%     end
%     obj.close
% setTag(t,'Compression',Tiff.Compression.None);
% setTag(t,'BitsPerSample',8);
% 
% setTag(t,'SampleFormat',Tiff.SampleFormat.UInt);
% setTag(t,'ExtraSamples',Tiff.ExtraSamples.Unspecified);

% setTag(t,'TileLength',32);
% setTag(t,'TileWidth',32);
% 

end

end








