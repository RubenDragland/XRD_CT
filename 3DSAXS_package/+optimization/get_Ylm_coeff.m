function Ylm_out = get_Ylm_coeff(ll,mm)
% Ylm = get_Ylm_coeff(l,m) returns precomputed coefficients for the spherical 
% harmonics version of polar order l and azimuthal order m.
%
% Note that for a given (l,m) if the orientation is optimized one needs
% also the coefficients for (l,m+1)
%
% Note that the order of rows has to be changed according to the order of
% the input coefficients, so here we change the order of rows (first index)
% to reflect the orders given in the input l and m. However, note that we
% dont change the order of the second coefficient, the column order. This
% is because the column order is linked in SAXS_tomo_3D_err_metric to the
% order of the power of cosines that we use in block_cos_theta_powers,
% this means that if we dont change the order in block_cos_theta_powers we
% also dont need to change here the order of the columns. For
% generalization to more complex spherical harmonics we have to preserve
% the fact that for each column here there should be a corresponding
% function of powers of sin(theta) and cos(theta) in block_cos_theta_powers
% in SAXS_tomo_3D_err_metric. How to do this when many coefficients are
% possible for l = 6 and many m's remains something to think about.
%
%   Inputs
% l     vector with polar order of spherical harmonics
% m     vector with azimuthal order of spherical harmonics
%
%   Outputs
% Ylm   vector with coefficients
%
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

if ~all(size(mm)==size(ll))
    error('l and m must be of the same size')
end

Ylm_out = zeros(numel(mm));

for ii = 1:numel(mm)
    m = mm(ii);
    l = ll(ii);
    
    if l == 0 && m == 0
        Ylm_out(ii,1) = 1/(2.*sqrt(pi));
    elseif l == 0 && m == 1
        Ylm_out(ii,1) = 0;
        
    elseif l == 2 && m == 0
        Ylm_out(ii,1:2) = sqrt(5/pi)*[-1/4, 3/4 ];
    elseif l == 2 && m == 1
        Ylm_out(ii,1:2) = sqrt(15/(2*pi))*[0, -1/2];
        
    elseif l == 4 && m == 0
        Ylm_out(ii,1:3) =  sqrt(1/pi)*[9/16, -45/8, 105/16];%3/16.*sqrt(1/pi).*(35.*(cos(theta)).^4 - 30.*(cos(theta)).^2 +3);
    elseif l == 4 && m == 1
        Ylm_out(ii,1:3) =  sqrt(5/pi)*[0, 9/8, -21/8];
        
    elseif l == 6 && m == 0
        Ylm_out(ii,1:4) =  sqrt(13/pi)*[-5/32, 105/32, -315/32, 231/32];%3/16.*sqrt(1/pi).*(35.*(cos(theta)).^4 - 30.*(cos(theta)).^2 +3);
    elseif l == 6 && m == 1
        Ylm_out(ii,1:4) =  sqrt(273/(2*pi))*[0, -5/16, 15/8, -33/16];
  
%     elseif l == 2 && m == 1
%         Ylm_out = -1/2.*sqrt(15/(2.*pi)).*sin(theta).*cos(theta).*exp(1i.*phi);
%         
%     elseif l == 4 && m == 1
%         Ylm_out =  -3/8.*sqrt(5/pi).*exp(1i.*phi).*sin(theta).*(7.*(cos(theta)).^3-3.*cos(theta));
%         
    else
        error('Combination of (l,m) = (%d,%d) is not available as precomputed quantity')
    end
end

return