# Repository for the small-angle scattering tensor tomography (SASTT)

## Software requirements
A version of matlab starting from **2016b**
Additional dependences:
- Parallel Toolbox
- Signal Processing Toolbox (window3.m)
- Image Processing Toolbox (step4_airtransmission_2Dmask.m)

## Tips on running scripts
By default, a parallel pool (if none exists yet) starts with a number of workers equal to the number of cores on the machine, **ignoring** "Preferred number of workers in a parallel pool" preference. To use a specific number of cores, initiate parallel pool **before** running scripts with
```matlab
parpool('local', desired_number_of_cores);
```

## References
* M. Liebi, M. Georgiadis, A. Menzel, P. Schneider, J. Kohlbrecher, O. Bunk, and M. Guizar-Sicairos, “Nanostructure surveys of macroscopic specimens by small-angle scattering tensor tomography,” *Nature* **527**, 349-352 (2015). (doi:10.1038/nature16056)
* M. Liebi, M. Georgiadis, J. Kohlbrecher, M. Holler, J. Raabe, I. Usov, A. Menzel, P. Schneider, O. Bunk and M. Guizar-Sicairos, "Small-angle X-ray scattering tensor tomography: model of the three-dimensional reciprocal-space map, reconstruction algorithm and angular sampling requirements," *Acta Cryst.* **A74**, 12-24 (2018). (doi:10.1107/S205327331701614X)
 
## Latest license
```
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
```