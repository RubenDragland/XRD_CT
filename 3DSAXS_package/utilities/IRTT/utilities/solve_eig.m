%% vecs = solve_eig(tensor, which_eig, thres_low,thres_high)
% Calculate eigen vectors per voxel from tensor model
%
% Inputs:
%   tomotensor  reconstructed tensor model
%   which_eig   which eigenvertor to solve
%
% Above is required for user input.
%   thres_low   low threshold for vector length
%   thres_high  high threshold for vector length

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

function vecs = solve_eig(tensor, which_eig, thres_low,thres_high)
    if nargin < 4
        thres_high = inf;
        if nargin < 3
            thres_low = 0;
            if nargin < 2
                which_eig = 1;
            end
        end
    end
    which_eig = 4 - which_eig;
    nx = size(tensor,2);
    ny = size(tensor,1);
    nz = size(tensor,3);
    vecs = zeros(ny,nx,nz,3);
    for i = 1:ny
        for j = 1:nx
            for k = 1:nz
                xx = tensor(i,j,k,1); 
                yy = tensor(i,j,k,2); 
                zz = tensor(i,j,k,3); 
                xy = tensor(i,j,k,4); 
                xz = tensor(i,j,k,5); 
                yz = tensor(i,j,k,6);
                T = [xx,xy,xz;xy,yy,yz;xz,yz,zz];
                [eigvecs, eigvals] = eig(T);
                eigvals = abs(eigvals);
                
                if (eigvals(which_eig,which_eig) > thres_low)&&(eigvals(which_eig,which_eig) < thres_high)
                    vecs(i,j,k,:) = eigvecs(:,which_eig)*eigvals(which_eig, which_eig);
                end
            end
        end
    end
end
