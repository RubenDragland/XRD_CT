% function projection = filter_SASTT_diode(projection,n,filtertype)
%
% Receives SASTT projections and denoises the transmission correction. This is 
% achieved by multiplying by the diode transmission (so undoing the transmission
% compensation), filtering the diode transmissiona and finally dividing the 
% filtered transmission
%  Inputs
% n             Size of the kernel
% filtertype    Either 'mean' or 'median'


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

function projection = filter_SASTT_diode(projection,n,filtertype)

for ii = 1:numel(projection)
    
    trans = projection(ii).diode;
    
    trans = padarray(trans,[n n],'replicate','both');
        
    switch lower(filtertype)
        case 'median'
            trans = medfilt2(trans,[n n],'symmetric');
        case 'mean'
            kernel = ones(n);
            kernel = kernel/sum(kernel(:));
            trans = conv2(trans,kernel,'same');
    end
    
    trans = trans(n+1:end-n,n+1:end-n);
    
    for jj = 1:size(projection(ii).data,3)
        projection(ii).data(:,:,jj) = projection(ii).data(:,:,jj).*projection(ii).diode./trans;
    end
end