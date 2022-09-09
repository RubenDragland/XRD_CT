function [ plothandle ] = plot_spherical_harmonic_squared(orientation,l,m,a,plotradius)
% [ plothandle ] = plot_spherical_harmonic(orientation,l,m,a,plotradius)
% Plots an addition of weighted spherical harmonics with an arbitrary axis
% orientation
% Inputs
%   orientation = [theta phi], Rotation of axes is performed first by theta
%               along y and then by phi along z.
%   l   Vector of polar orders
%   m   Vector of azimuthal orders
%   a   Vector of coefficients
%   plotradius     If false the function is shown on a sphere, if true the
%                   radius is proportional to the function
% Outputs
%   plothandle  (optional output)
%
%
%   Note that orientation defines a gendataeral orientation of z, but not for
%   (x,y), if the azimuthal orientation of this function is important the
%   rotation matrices below must be generated in the same way in the
%   optimization
%
%   Manuel Guizar 19 Nov 2013
%
%   adjust notations
%
%   Marianne Liebi 23 May 2014
%
%   plot  Ymn_squared = abs((Ymn).^2) which is the result of the
%   optimization
    
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

theta_struct = orientation(1);
phi_struct = orientation(2);

% % Rot_str is the combination of 1. Rotation around y with theta_struct
    % and 2. Rotation around z with theta_struct
    % this rotation is needed from the spherical harmonics coordinate system (main
    % orientation axis along z-axis to the object coordinate system
    Rot_y_theta_struct = [cos(theta_struct) 0 -sin(theta_struct) ; 0 1 0 ; sin(theta_struct) 0 cos(theta_struct)];
    Rot_z_phi_struct = [cos(phi_struct) sin(phi_struct) 0 ; -sin(phi_struct) cos(phi_struct) 0 ; 0 0 1];
    Rot_str = Rot_y_theta_struct * Rot_z_phi_struct;
    Rot = Rot_str; %here the rotation is only from spherical harmonics to coordinate system, not in beamline coordi
    
% Assingment of convenient variables
alphai = Rot(1,1);
alphaj = Rot(1,2);
alphak = Rot(1,3);
betai  = Rot(2,1);
betaj  = Rot(2,2);
betak  = Rot(2,3);
gammai = Rot(3,1);
gammaj = Rot(3,2);
gammak = Rot(3,3);


%  Create a grid of angles
delta = pi/40;
theta = 0 : delta : pi; % altitude
phi = 0 : 2*delta : 2*pi; % azimuth
[phi,theta] = meshgrid(phi,theta);

x_sh = (   alphai*sin(theta).*cos(phi) + alphaj*sin(theta).*sin(phi) + alphak*cos(theta)   );
y_sh = (   betai*sin(theta).*cos(phi) + betaj*sin(theta).*sin(phi) + betak*cos(theta)   );
z_sh = (   gammai*sin(theta).*cos(phi) + gammaj*sin(theta).*sin(phi) + gammak*cos(theta)   );

theta_sh = acos(z_sh);
phi_sh = atan2(y_sh,x_sh); 

Ymn = zeros(size(theta_sh));
for jj = 1:numel(l)
    Ymn = Ymn + a(jj)*optimization.spherical_harm(l(jj),m(jj),theta_sh,phi_sh);
end

%Ymn = real(Ymn); % The real part of the imaginary azimuthal component is physical
 Ymn_squared = abs((Ymn).^2); %in order to avoid complex and negative intensities
    %     data_synt(:,ii) = real(Ymn_2D);

switch plotradius
    case 0
        X = sin(theta).*cos(phi);
        Y = sin(theta).*sin(phi);
        Z = cos(theta);
    case 1
        X =  Ymn_squared.*sin(theta).*cos(phi);
        Y =  Ymn_squared.*sin(theta).*sin(phi);
        Z =  Ymn_squared.*cos(theta);
end

plothandle =  surf(X,Y,Z, Ymn_squared);
% % surf(X,Y,Z,Theta)
% % surf(X,Y,Z,Phi)
axis equal
colorbar
xlabel('x')
ylabel('y')
zlabel('z')
% % view(-90,0)


end

