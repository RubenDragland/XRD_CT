% w = window3(N,M,L,w_type)
%
% Creates a 3D lowpass filtering window using the window function in Matlab, it takes
% as arguments scalar dimensions of the array (N,M,L) and the
% filter type
% ---------------------------------------------------------------------
%     w_type can be: 
%     @bartlett       - Bartlett window.
%     @barthannwin    - Modified Bartlett-Hanning window. 
%     @blackman       - Blackman window.
%     @blackmanharris - Minimum 4-term Blackman-Harris window.
%     @bohmanwin      - Bohman window.
%     @chebwin        - Chebyshev window.
%     @flattopwin     - Flat Top window.
%     @gausswin       - Gaussian window.
%     @hamming        - Hamming window.
%     @hann           - Hann window.
%     @kaiser         - Kaiser window.
%     @nuttallwin     - Nuttall defined minimum 4-term Blackman-Harris window.
%     @parzenwin      - Parzen (de la Valle-Poussin) window.
%     @rectwin        - Rectangular window.
%     @taylorwin      - Taylor window.
%     @tukeywin       - Tukey window.
%     @triang         - Triangular window.
%
%   Example: 
%   To compute windowed 2D fFT
%   [r,c,l]=size(img);
%   w=window3(r,c,l,@hamming);
% 	fftn(img.*w);
%
% Marianne Liebi, Oct 26, 2016
%
% This function was extended from the original submission to Mathworks by
% Disi A ---- May,16, 2013
% http://ch.mathworks.com/matlabcentral/fileexchange/43827-two-dimensional-window-generator--2d-window-

% Copyright (c) 2013, Disi A
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% Copyright (c) 2016, Paul Scherrer Institut, Marianne Liebi
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

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


function w=window3(N,M,L,w_func)

wc=window(w_func,N);
wr=window(w_func,M);
wf=window(w_func,L);
[maskr,maskc, maskf]=meshgrid(wr,wc,wf);

%maskc=repmat(wc,1,M); Old version
%maskr=repmat(wr',N,1);

w=maskr.*maskc.*maskf;

end