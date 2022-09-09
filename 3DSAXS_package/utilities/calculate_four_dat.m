% function [four_dat] = calculate_four_dat(four_dat, proj)
%
% Calculates parameters needed for plotting

%*-------------------------------------------------------------------------------------*
%|                                                                                     |
%|  Except where otherwise noted, this work is licensed under a                        |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0                          |
%|  International (CC BY-NC-SA 4.0) license.                                           |
%|                                                                                     |
%|  Copyright (c) 2018 by Paul Scherrer Institute (http://www.psi.ch)                  |
%|                                                                                     |
%|      Author: CXS group, PSI                                                         |
%*------------------------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If this code, or subfunctions or parts of it, is used for research in a
%   publication or if it is fully or partially rewritten for another
%   computing language the authors and institution should be acknowledged
%   in written form and additionally you should cite:
%     M. Liebi, M. Georgiadis, A. Menzel, P. Schneider, J. Kohlbrecher,
%     O. Bunk, and M. Guizar-Sicairos, “Nanostructure surveys of
%     macroscopic specimens by small-angle scattering tensor tomography,”
%     Nature 527, 349-352 (2015).   (doi:10.1038/nature16056)
% and
%     M. Liebi, M. Georgiadis, J. Kohlbrecher, M. Holler, J. Raabe, I. 
%     Usov, A. Menzel, P. Schneider, O. Bunk and M. Guizar-Sicairos,
%     "Small-angle X-ray scattering tensor tomography: model of the 
%     three-dimensional reciprocal-space map, reconstruction algorithm 
%     and angular sampling requirements," Acta Cryst. A74, 12-24 (2018). 
%     (doi:10.1107/S205327331701614X)
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its
%    proper use and the correctness of the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [four_dat] = calculate_four_dat(four_dat, proj)

% calculate fourier coeffiecients
f1_amp = [];
f2_amp = [];
f2_phase = [];
for ii = 1:size(proj, 1)
    for jj = 1:size(proj, 2)
        point = squeeze(proj(ii, jj, :));
        [f1_amp1, f2_amp1, f2_phase1] = fourier_coefficients(point);
        f1_amp(ii, jj, :) = f1_amp1;
        f2_amp(ii, jj, :) = f2_amp1;
        f2_phase(ii, jj, :) = f2_phase1;
    end
end
% Convert the phase of the fft coefficient to orientation in degrees
% flag for ignoring the sign of f2_phase
ignore_f2p_sign = 0;
phi_det = four_dat.integ.phi_det;
orientation = phase2orientation(f2_phase,phi_det,ignore_f2p_sign);
degree_orient = f2_amp ./ (f1_amp + 1e-6); % to avoid a division by zero

% replace data in four_dat
four_dat.f1_amp = f1_amp;
four_dat.f2_amp = f2_amp;
four_dat.f2_phase = f2_phase;
four_dat.orientation = orientation;
four_dat.degree_orient = degree_orient;
        