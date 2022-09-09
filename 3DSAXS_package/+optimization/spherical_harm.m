function Ylm = spherical_harm(l,m,theta,phi)
% Ylm = spherical_harm(l,m,theta,phi) computes the orthonormal version of
% the spherical harmonic of polar order l and azimuthal order m in polar
% coordinate theta and azimuthal coordinate phi.
%
% Theta and phi can be at most a three dimensional array and they must have
% the same dimensions.
%
% This implementation uses native matlab legendre function that computes
% all m = [0 l] orders. If you are looping over m outside this code then
% you are waisting precious computation
%
% Manuel Guizar - Nov 19 2013
%
% write Ylm explicitly for l=0, m=0; l=2, m=0 and l=4, m=0
% values have been checked
%
% Marianne Liebi 23 May 2014
%
% include l=0, m=0+1 to run the loop in the gradient without crashing
% the following check is therefor commented:
% if abs(m)>l
%     error('Azimuthal order m cannot be larger than l')
% end
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

if ~all(size(theta) == size(phi))
    error('theta and phi must have the same dimensions')
end



if l == 0 && ...
        m == 0
    Ylm = 1/(2.*sqrt(pi));
    
elseif l == 0 && ...
        m == 1
    Ylm = 0; %%%this is included to run the loop without crashing for l=0 and m=0+1
    
elseif l == 2 && ...
        m == 0
    Ylm = 1/4.*sqrt(5/pi).*(3.*(cos(theta)).^2 -1);
    
elseif l == 4 && ...
        m == 0
    Ylm = 3/16.*sqrt(1/pi).*(35.*(cos(theta)).^4 - 30.*(cos(theta)).^2 +3);
    
elseif l == 2 &&...
    m == 1
    Ylm = -1/2.*sqrt(15/(2.*pi)).*sin(theta).*cos(theta).*exp(1i.*phi);
    
elseif l == 4 &&...
    m == 1
    Ylm = -3/8.*sqrt(5/pi).*exp(1i.*phi).*sin(theta).*(7.*(cos(theta)).^3-3.*cos(theta));
    
    
else
    Ylm = legendre(l,cos(theta));
    % Because of matlab implementation now I have in the first index all orders
    % from m = 0 to l
    % This could be exploited later better for looping over m,
    % for now I just choose one that I need
    if l~=0
        switch numel(size(theta))
            case 1
                Ylm = squeeze(Ylm(abs(m)+1,:));
            case 2
                Ylm = squeeze(Ylm(abs(m)+1,:,:));
            case 3
                Ylm = squeeze(Ylm(abs(m)+1,:,:,:));
            otherwise
                error('Arrays theta and phi have more than 3 dimensions')
        end
    end
    % Checking for negative m
    if m < 0
        Ylm = (-1)^m*factorial(l+m)*Ylm/factorial(l-m);
        %     Ylm = (-1)^m*Ylm;
    end
    
    
    % Normalization
    Ylm = sqrt(   (2*l+1)*factorial(l-m)/((4*pi)*(factorial(l+m)))   )*Ylm.*exp(1i*m*phi);
end

return