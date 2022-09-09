%% tcorr=proj_tensor_correct(tomotensor,R,diff,projsize)
% Correct the tensor model according to the differential map.
%
% Inputs:
%   tomotensor  reconstructed tensor model to be corrected
%   R           spatial rotation matrix
%   diff        differential map between measured and simulated projections
%   projsize    size of the projection generated

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

function tcorr=proj_tensor_correct(tomotensor,R,diff,projsize)

    modelsize = size(tomotensor);
    if nargin < 4
        projsize = [modelsize(4), modelsize(1), modelsize(2)];
    end
    xgrid = (-modelsize(1)/2+1 ):(modelsize(1)/2);
    ygrid = (-modelsize(2)/2+1 ):(modelsize(2)/2);
    zgrid =  (-modelsize(3)/2+1) : (modelsize(3)/2);
    [Y,X,Z] = meshgrid(ygrid, xgrid, zgrid);
    
    Xp = R(2,2)*X + R(2,1)*Y + R(2,3)*Z + projsize(2)/2;
    Yp = R(1,2)*X + R(1,1)*Y + R(1,3)*Z + projsize(3)/2;
    
    length = min(modelsize(2)/abs(R(2,2)),modelsize(1)/abs(R(3,2)));

    tcorr = projection_tensor_correct(tomotensor,modelsize,projsize,Xp,Yp,diff/length); 
    
    
    
    
   
%     tomo = tomotensor;
%     modelsize = size(tomo);
%     radius = (modelsize(1)-1)/2;
%     if nargin<6
%         projsize = [modelsize(3),modelsize(2)];
%     end
% 
%     rotx = rotx/180*pi;
%     roty = roty/180*pi;
%     
%     stepx = cos(rotx)*cos(roty);
%     stepy = -cos(rotx)*sin(roty);
%     stepz = sin(rotx);
%     
%     steppixx = sin(roty);
%     steppixy = cos(roty);
%     
%     steppixz = 1/cos(rotx);
%     
%     startx0 = radius*(1-cos(roty)-sin(roty));
%     starty0 = radius*(1+sin(roty)-cos(roty));
%     startz = (modelsize(3)-projsize(1)/cos(rotx))/2-radius*tan(rotx);

end
