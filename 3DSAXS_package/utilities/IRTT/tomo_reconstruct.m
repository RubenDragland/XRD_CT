%% [tomogram,projection]=tomo_reconstruct(filename,ifshow)
% Tomographic reconstruction with filtered back projection per slice using
% only projections from 0 tilt angle.
% Also loads data from saved mat file.
%
% Inputs:
%   filename    name of the saved mat file, or the projection structure
%   containing SAXS data
%   ifshow      if to show the reconstructed slices, =1 by default

% Copyright 2017 Zirui Gao
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [tomogram,projection]=tomo_reconstruct(filename,ifshow)
    
    if ischar(filename)
        projection = load_projections(filename);
    else
        projection = filename;
    end
    
    if nargin<2
        ifshow=1;
    end
    
    N = length(projection);
    
    % modify until here.
    i=1;
    while (i<=N)&&(projection(i).rot_x==0)
        roty(i)=projection(i).rot_y;
        projs(i,:,:)=mean(projection(i).data,3);
        i=i+1;
    end
    tomogram = zeros(size(projs,2),size(projs,3),size(projs,3));
    for i = 1:size(projs,2)
        sinotrans = squeeze(projs(:,i,:))';
        sinotrans = utils.shiftwrapbilinear(sinotrans,0,0);
        tomogram(i,:,:)  = iradon(sinotrans,-roty,'linear','Hamming',size(projs,3),1)';
        if ifshow==1
           figure(114);
           imagesc(squeeze(tomogram(i,:,:)))
           axis xy equal tight;    colormap bone;    colorbar
           title(['Slice ' num2str(i)])
           pause(0.2)
        end
    end
end