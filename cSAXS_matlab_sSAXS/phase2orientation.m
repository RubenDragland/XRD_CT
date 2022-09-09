% function [orientation] = phase2orientation(f2_phase,phi_det,ignore_f2p_sign)
%
% Converts the phase of Fourier coefficients to orientation in degrees.
% Please note this is very particular for the case of sSAXS, because the
% Fourier coefficient covers a phase of 2*pi radians for 180 degrees
% orientation, basically assuming that the pattern is centrosymmetrical.
%
% orientation = -phase/2 + phi_det(1)/2    (1)
%
% Inputs
% f2_phase  Phase of the complex FFT coefficient
% phi_det   Angular center of the segments, just uses the first element and
%           it assumes that the angular spacing is equal
% ignore_f2p_sign   Some debugging flag that ignores the sign of the
%                   orientation

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%*-----------------------------------------------------------------------*
%
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing la this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a
%   publication or if it is fully or partially rewritten for another
%   computing language the authors and institution should be acknowledged
%   in written form in the publication: “Data processing was carried out
%   using the “cSAXS scanning SAXS package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.”
%   Variations on the latter text can be incorporated upon discussion with
%   the CXS group if needed to more specifically reflect the use of the package
%   for the published work.
%
% Additionally, any publication using the package, or any translation of the
%     code into another computing language should cite:
%    O. Bunk, M. Bech, T. H. Jensen, R. Feidenhans'l, T. Binderup, A. Menzel
%    and F Pfeiffer, “Multimodal x-ray scatter imaging,” New J. Phys. 11,
%    123016 (2009). (doi: 10.1088/1367-2630/11/12/123016)
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%
% This code and subroutines are part of a continuous development, they
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its
%    proper use and the correctness of the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [orientation] = phase2orientation(f2_phase,phi_det,ignore_f2p_sign)
% offset the zero-point of the rotation from the start of the first 
% segment to its center
% four_dat.f2_phase = four_dat.f2_phase - pi/8;
f2_phase = f2_phase - phi_det(1)*pi/90; % 90 instead of 180 because of Eq. (1) above
if (ignore_f2p_sign)
    f2_phase = abs(f2_phase);
end

% phase from -90 to 90 degree
orientation = -f2_phase / pi * 90.0; % 90 instead of 180 because of Eq. (1) above
% Wrap orientation from -90 to 90
ind_wrap = find(orientation < -90.0);
orientation(ind_wrap) = orientation(ind_wrap) + 180.0;
ind_wrap = find(orientation >= 90.0);
orientation(ind_wrap) = orientation(ind_wrap) - 180.0;

end