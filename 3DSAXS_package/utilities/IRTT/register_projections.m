%% projection_aligned=register_projections(projection, tomogram, usemask,checkonly)
% Register all the projections to the model from tomographic reconstruction
% DO NOT run this script again if the alignment is already done with ASTRA.
%
% Inputs:
%   projection  structure containing data of SAXS projections.
%   tomogram    reconstructed scalar tomogram of the sample
%
% Above is required for user input.
%   usemask     if to mask the projections with projection.window_mask
%   before alignment, =0 by default
%   checkonly   only check the alignment results with out doing alignment,
%   =0 by default

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


function projection_aligned=register_projections(projection, tomogram, usemask,checkonly)

    if nargin<3
        usemask=0;
    end
    if nargin<4
        checkonly=0;
    end
    
    numangles=size(projection,2);
%     delta_proj=zeros(numangles,2);
    for i=1:numangles
        praw=mean(projection(i).data,3);
        psim=proj_simul(tomogram,projection(i).Rot_exp,size(praw));
        
        if usemask
            praw=praw.*projection(i).window_mask;
        end
        if checkonly
            delta=[projection(i).dy projection(i).dx];
            figure(115); imagesc(psim); axis equal tight; ca=caxis; title(sprintf('rotx %d roty %.1f',projection(i).rot_x,projection(i).rot_y));
            figure(116); imagesc(utils.shiftwrapbilinear(praw,-delta(1),-delta(2))); axis equal tight; caxis(ca); title(sprintf('rotx %d roty %.1f',projection(i).rot_x,projection(i).rot_y));
%             disp(delta);
            pause(0.1);
        else
            try
                [~,~,delta]=utils.registersubimages_2(praw,psim,[],[],[],[],100,0,1);
            catch err
                disp ('Error in registersubimages_2!');
                delta=[0 0];
            end
            if 1
                %             praw=proj_translation(praw,delta);
                figure(115); imagesc(psim); axis equal tight; ca=caxis; title(sprintf('rotx %d roty %.1f',projection(i).rot_x,projection(i).rot_y));
                figure(116); imagesc(utils.shiftwrapbilinear(praw,-delta(1),-delta(2))); axis equal tight; caxis(ca); title(sprintf('rotx %d roty %.1f',projection(i).rot_x,projection(i).rot_y));
                %             disp(delta);
                pause(0.1);
            end
            projection(i).dy=delta(1);
            projection(i).dx=delta(2);
%         delta_proj(i,:)=delta(:);
        end
    end
    projection_aligned = projection;
%       save([path_recon '/delta_proj.mat'],'delta_proj');
end


% set(gca,'position',[0 0 1 1],'units','normalized');