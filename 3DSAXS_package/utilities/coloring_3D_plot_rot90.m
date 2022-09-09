function coloring_3D_plot_rot90(theta_struct_slice,phi_struct_slice,deg_orient,fig_new)
% Create spheres as colormaps

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


% Make a matrix for 3D coloring
clear colormatrix
colormatrix(1,:,:) = fftshift([colormap(hsv(81));colormap(hsv(81))],1);
colormatrix = repmat(colormatrix,[80 1 1]);
colormatrix = rgb2hsv(colormatrix);
colormatrix(:,:,2) = repmat([linspace(0,1,40).';linspace(1,0,40).'],[1 162]);
colormatrix = hsv2rgb(colormatrix);
[colormatrix,map] = rgb2ind(colormatrix,2000);
colormatrix = im2double(colormatrix,'indexed');


% figure(5)
% imshow(colormatrix,map)
% 
% figure(4)
% [x,y,z] = sphere;
% % surf(x,y,z)  % sphere centered at origin
% surf(x,y,z,colormatrix,...
%    'FaceColor','texturemap',...  
%    'EdgeColor','none',...
%    'CDataMapping','direct')
% colormap(map)
% daspect([1 1 1])
% view(2)
% grid off

figure(6)
[x,y,z] = sphere;
% surf(x,y,z)  % sphere centered at origin
surf(x,y,z,colormatrix,...
   'FaceColor','texturemap',...  %'EdgeColor','none',...
   'CDataMapping','direct')
colormap(map)
daspect([1 1 1])
view(0,75)
grid off
% function which converts 3D orientation in color,mapped on one hemisphere (direction of vector irrelevant) 
%hue:phi  
%saturation: theta
%value: optional, degree of orientation
%
%input: values are nxm matrices, all of the same size
%
%ouput: plot with colors representing the orientation
%the polar angles points towards qou (z). If the orientation vector points
%in the xy plane then you see bright colors, as the vector points towards
%you (positive z) then the color. the colors going twice around the
%azimuthal to avoid discontinuities, but same color can correspond to two
%different orientations





unitx=sin(theta_struct_slice).*cos(phi_struct_slice);
unity=sin(theta_struct_slice).*sin(phi_struct_slice);
unitz=cos(theta_struct_slice);

    % All pointing at the positive z hemisphere
    unitx = unitx.*sign(unitz+eps);
    unity = unity.*sign(unitz+eps);
    unitz = unitz.*sign(unitz+eps);
    % Conversion to traditional polar coordinates with z as main axis
    phiconv   = atan2d(unity,unitx); % In degrees [-180 180]
    
   
        phiconv(phiconv<0) = phiconv(phiconv<0)+180; % Ensure only positive angles within 0 and 180
 
        phiconv(phiconv>180) = phiconv(phiconv>180)-180;
 
    thetaconv = atan2d(sqrt(unitx.^2+unity.^2),unitz);
    % Put results on individual arrays for easy inspection
%     Phiconv(coinc_map1_improved{2}(ii,1)+1, coinc_map1_improved{2}(ii,2)+1) = phiconv;
%     Thetaconv(coinc_map1_improved{2}(ii,1)+1, coinc_map1_improved{2}(ii,2)+1) = thetaconv;
    % Hue to azimuthal angle
    phiconv = phiconv/180;
    angles3D_hsv(:,:,1) = phiconv;
    % Saturation coded to theta orientation, white points towards you
    thetaconv = thetaconv/90; % When theta = 0 saturation is zero (white when it points at you)
    angles3D_hsv(:,:,2) = thetaconv;
    % Full intensity on coincident points
    angles3D_hsv(:,:,3) = deg_orient;


angles3D = hsv2rgb(angles3D_hsv);

% figure(1); 
% % plot(coinc_map1_improved{2}(:,2),coinc_map1_improved{2}(:,1),'.b')
% % imagesc(angles3D_hsv(:,:,2))
% % imagesc(phi)
% % imagesc(theta)
% % imagesc(Phiconv)
% imagesc(thetaconv)
% axis equal tight
% colormap bone
% colorbar
if fig_new
figure;
image(rot90(angles3D,3))
axis xy equal tight
else
figure(100);
image(rot90(angles3D,3))
axis xy equal tight   
end  



end

