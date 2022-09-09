function outerror =  errorplot(argin)
% Build an errormetric plot vs iteration number
% e = errorplot         Clears the persistent variable
% e = errorplot(x)      Appends x to the persistent variable
% e = errorplot([])     Only reads the persistent variable

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

persistent errormetric,

if exist('argin') == 0,
    errormetric = [];
else
    errormetric = [errormetric ; argin];
end

outerror = errormetric;