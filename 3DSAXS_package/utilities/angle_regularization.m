function [Reg, grad_theta_reg, grad_phi_reg] = angle_regularization(...
    theta_tomo, phi_tomo, sin_theta_struct, cos_theta_struct, sin_phi_struct, cos_phi_struct, ...
    ny, nx, nz)
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

if nargout ~= 1 && nargout ~= 3
    error('angle_regularization:nargout_check', ...
        'Function supports calls with only 1 or 3 output arguments');
end


sin_theta_tomo = reshape(sin_theta_struct, [ny,nx,nz]);
cos_theta_tomo = reshape(cos_theta_struct, [ny,nx,nz]);
sin_phi_tomo = reshape(sin_phi_struct,     [ny,nx,nz]);
cos_phi_tomo = reshape(cos_phi_struct,     [ny,nx,nz]);

dot_ind1   = ones([ny,nx,nz]);
dot_ind2   = ones([ny,nx,nz]);
dot_ind3   = ones([ny,nx,nz]);
one_matrix = ones([ny,nx,nz]);

dot_ind1(1:end-1,:,:) = ...
    sin_theta_tomo(1:end-1,:,:).*sin_theta_tomo(2:end,:,:).*cos(phi_tomo(1:end-1,:,:)-phi_tomo(2:end,:,:)) + ...
    cos_theta_tomo(1:end-1,:,:).*cos_theta_tomo(2:end,:,:);

dot_ind2(:,1:end-1,:) = ...
    sin_theta_tomo(:,1:end-1,:).*sin_theta_tomo(:,2:end,:).*cos(phi_tomo(:,1:end-1,:)-phi_tomo(:,2:end,:)) + ...
    cos_theta_tomo(:,1:end-1,:).*cos_theta_tomo(:,2:end,:);

dot_ind3(:,:,1:end-1) = ...
    sin_theta_tomo(:,:,1:end-1).*sin_theta_tomo(:,:,2:end).*cos(phi_tomo(:,:,1:end-1)-phi_tomo(:,:,2:end)) + ...
    cos_theta_tomo(:,:,1:end-1).*cos_theta_tomo(:,:,2:end);

Reg = sum(sum(sum(one_matrix-abs(dot_ind1)+one_matrix-abs(dot_ind2)+one_matrix-abs(dot_ind3))));

if nargout == 3
    grad_theta_reg = zeros([ny,nx,nz]);
    grad_phi_reg   = zeros([ny,nx,nz]);
    
    % Define utility variables
    sin_theta_cent = sin_theta_tomo(2:end-1,2:end-1,2:end-1);
    cos_theta_cent = cos_theta_tomo(2:end-1,2:end-1,2:end-1);
    sin_phi_cent = sin_phi_tomo(2:end-1,2:end-1,2:end-1);
    cos_phi_cent = cos_phi_tomo(2:end-1,2:end-1,2:end-1);
    
    sin_theta_up = sin_theta_tomo(1:end-2,2:end-1,2:end-1);
    sin_theta_low = sin_theta_tomo(3:end,2:end-1,2:end-1);
    sin_theta_left = sin_theta_tomo(2:end-1,1:end-2,2:end-1);
    sin_theta_right = sin_theta_tomo(2:end-1,3:end,2:end-1);
    sin_theta_front = sin_theta_tomo(2:end-1,2:end-1,1:end-2);
    sin_theta_back = sin_theta_tomo(2:end-1,2:end-1,3:end);
    
    cos_theta_up = cos_theta_tomo(1:end-2,2:end-1,2:end-1);
    cos_theta_low = cos_theta_tomo(3:end,2:end-1,2:end-1);
    cos_theta_left = cos_theta_tomo(2:end-1,1:end-2,2:end-1);
    cos_theta_right = cos_theta_tomo(2:end-1,3:end,2:end-1);
    cos_theta_front = cos_theta_tomo(2:end-1,2:end-1,1:end-2);
    cos_theta_back = cos_theta_tomo(2:end-1,2:end-1,3:end);
    
    sin_phi_up = sin_phi_tomo(1:end-2,2:end-1,2:end-1);
    sin_phi_low = sin_phi_tomo(3:end,2:end-1,2:end-1);
    sin_phi_left = sin_phi_tomo(2:end-1,1:end-2,2:end-1);
    sin_phi_right = sin_phi_tomo(2:end-1,3:end,2:end-1);
    sin_phi_front = sin_phi_tomo(2:end-1,2:end-1,1:end-2);
    sin_phi_back = sin_phi_tomo(2:end-1,2:end-1,3:end);
    
    cos_phi_up = cos_phi_tomo(1:end-2,2:end-1,2:end-1);
    cos_phi_low = cos_phi_tomo(3:end,2:end-1,2:end-1);
    cos_phi_left = cos_phi_tomo(2:end-1,1:end-2,2:end-1);
    cos_phi_right = cos_phi_tomo(2:end-1,3:end,2:end-1);
    cos_phi_front = cos_phi_tomo(2:end-1,2:end-1,1:end-2);
    cos_phi_back = cos_phi_tomo(2:end-1,2:end-1,3:end);
    
    sign_dot_ind1 = sign(dot_ind1(2:end-1,2:end-1,2:end-1));
    sign_dot_ind2 = sign(dot_ind2(2:end-1,2:end-1,2:end-1));
    sign_dot_ind3 = sign(dot_ind3(2:end-1,2:end-1,2:end-1));
    
    exprA = sign_dot_ind1.*(cos_theta_up+cos_theta_low) + ...
        sign_dot_ind2.*(cos_theta_left+cos_theta_right) + ...
        sign_dot_ind3.*(cos_theta_front+cos_theta_back);
    
    exprB = sign_dot_ind1.*(sin_theta_up.*cos_phi_up + sin_theta_low.*cos_phi_low) + ...
        sign_dot_ind2.*(sin_theta_left.*cos_phi_left + sin_theta_right.*cos_phi_right) + ...
        sign_dot_ind3.*(sin_theta_front.*cos_phi_front + sin_theta_back.*cos_phi_back);
    
    exprC = sign_dot_ind1.*(sin_theta_up.*sin_phi_up + sin_theta_low.*sin_phi_low) + ...
        sign_dot_ind2.*(sin_theta_left.*sin_phi_left + sin_theta_right.*sin_phi_right) + ...
        sign_dot_ind3.*(sin_theta_front.*sin_phi_front + sin_theta_back.*sin_phi_back);
    
    %the gradient of the regularization term in respect to theta
    grad_theta_reg(2:end-1,2:end-1,2:end-1) = ...
        sin_theta_cent.*exprA - cos_theta_cent.*(cos_phi_cent.*exprB + sin_phi_cent.*exprC);
    
    %the gradient of the regularization term in respect to phi
    grad_phi_reg(2:end-1,2:end-1,2:end-1) = ...
        sin_theta_cent.*(sin_phi_cent.*exprB - cos_phi_cent.*exprC);
end
end

