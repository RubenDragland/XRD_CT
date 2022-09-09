%% B_segs = calc_B_segs(projection,scatter_direction)
% Calculate the B-vectors (name from diffusion MRI), which is the direction
% of structure probed in each segment of each projection. 
%
% Inputs:
%   projection  structure containing data of SAXS projections.
%   scatter_direction   =1 to assume scattering parallel to the
%   structural model (required for combination with SASTT). =0 by default.

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


function [B_segs] = calc_B_segs(projection,scatter_direction)


if nargin<2
    scatter_direction=0;
end

numsegments = size(projection(1).data, 3);
if numsegments == 8
    anglestep = pi/numsegments;
    segangles = [anglestep/2:anglestep:pi] +(1-scatter_direction)*pi/2;
elseif numsegments > 8
    anglestep = 2*pi/numsegments;
    segangles = [anglestep/2:anglestep:2*pi];
end

p(1:numsegments,1) = sin(segangles);
p(1:numsegments,2) = cos(segangles);
p(1:numsegments,3) = 0;
B_vecs = zeros(length(projection),numsegments,3);
B_segs = zeros(length(projection),numsegments,6);
bvec = zeros(3,1);

for i = 1:length(projection)
%     rotx = deg2rad(projection(i).rot_x);
%     roty = deg2rad(projection(i).rot_y);
    
    for j = 1:numsegments
        bvec(:) = squeeze(p(j,:));
%         bvec =  [cos(rotx), 0,  sin(rotx);...
%                     0,      1,      0;...
%                 -sin(rotx), 0,  cos(rotx)] * bvec; % rotation around X and Y axis
%         bvec =  [cos(roty), sin(roty),  0;...
%                 -sin(roty), cos(roty),  0;...
%                     0,          0,      1] * bvec;    
        R=projection(i).Rot_exp([2,1,3],[2,1,3]);
        bvec=R\bvec;
        B_vecs(i,j,:) = bvec;
        % Expand the vector to V*V', to later multiply with the tensor
        B_segs(i,j,:) = [bvec(1)^2, bvec(2)^2, bvec(3)^2, ...
            2*bvec(1)*bvec(2), 2*bvec(1)*bvec(3), 2*bvec(2)*bvec(3)]; 
    end
end
end