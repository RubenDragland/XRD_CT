first_scan_nr = 272;
Nx=11;
nr_lines = 11;
centerx = 224;
centery = 98;
dark_field_radius = 20;
roi = 128;
detector_number = 1; % <1-Pilatus 2M, 2-Pilatus 300k, 3-Pilatus 100k, 4-Eiger 500k>


[trans,dpcx,dpcy,df] = beamline.stxm_online(first_scan_nr, nr_lines , 'Nx',Nx,'ROIdim',roi,'CenX', centerx, 'CenY', centery, ...
    'DarkFieldR', dark_field_radius,'FilenameValidMask',[],'DirPerLine',0 , ...
    'DetectorNumber', detector_number);

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