%% projection = proj_simul_tensor(tomotensor,R,B,projsize)
% Generate simulated SAXS projection from a tensor model.
%
% Inputs:
%   tomotensor  reconstructed tensor model
%   R           spatial rotation matrix
%   B           B matrix of the relative orientation
%   projsize    size of the projection to be generated

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


function projection = proj_simul_tensor(tomotensor,R,B,projsize)
    tomo = tomotensor;
    modelsize = size(tomo);
    if nargin < 4
        projsize = [size(B,1),modelsize(1), modelsize(2)];
    end
    xgrid = (-modelsize(1)/2+1 ):(modelsize(1)/2);
    ygrid = (-modelsize(2)/2+1 ):(modelsize(2)/2);
    zgrid =  (-modelsize(3)/2+1) : (modelsize(3)/2);
    [Y,X,Z] = meshgrid(ygrid, xgrid, zgrid);
    
    Xp = R(2,2)*X + R(2,1)*Y + R(2,3)*Z + projsize(2)/2;
    Yp = R(1,2)*X + R(1,1)*Y + R(1,3)*Z + projsize(3)/2;

    projection = projection_sim_tensor(tomotensor,modelsize,projsize,Xp,Yp,B); 
% this part is rewritten in C for speed
%     projection=zeros(106,70);
%     for sz=1:106
%         startx=startx0;
%         starty=starty0;
%         for sy=1:70
%             x=startx;
%             y=starty;
%             z=startz;
%             t=1;
%             while ((sqrt((x-34.5)^2+(y-34.5)^2)>34.5)||(z<0)||(z>105))&&(t<100)
%                 x=x+stepx;
%                 y=y+stepy;
%                 z=z+stepz;
%                 t=t+1;
%             end
%             while (sqrt((x-34.5)^2+(y-34.5)^2)<=34.5)&&(z>=0)&&(z<=105)
%                 x2=x+1; y2=y+1; z2=z+1;
%                 xf=fix(x2); yf=fix(y2); zf=fix(z2);
%                 xinp00 = tomo(xf,yf,zf)*(xf+1-x2) + tomo(xf+1,yf,zf)*(x2-xf);
%                 xinp10 = tomo(xf,yf+1,zf)*(xf+1-x2) + tomo(xf+1,yf+1,zf)*(x2-xf);
%                 xinp11 = tomo(xf,yf+1,zf+1)*(xf+1-x2) + tomo(xf+1,yf+1,zf+1)*(x2-xf);
%                 xinp01 = tomo(xf,yf,zf+1)*(xf+1-x2) + tomo(xf+1,yf,zf+1)*(x2-xf);
%                 yinp0 = xinp00*(yf-y2+1)+xinp10*(y2-yf);
%                 yinp1 = xinp01*(yf-y2+1)+xinp11*(y2-yf);
%                 out = yinp1*(z2-zf)+yinp0*(zf-z2+1);
% 
%                 xinp00 = getv(xf,yf,zf)*(xf+1-x2) + getv(xf+1,yf,zf)*(x2-xf);
%                 xinp10 = getv(xf,yf+1,zf)*(xf+1-x2) + getv(xf+1,yf+1,zf)*(x2-xf);
%                 xinp11 = getv(xf,yf+1,zf+1)*(xf+1-x2) + getv(xf+1,yf+1,zf+1)*(x2-xf);
%                 xinp01 = getv(xf,yf,zf+1)*(xf+1-x2) + getv(xf+1,yf,zf+1)*(x2-xf);
%                 yinp0 = xinp00*(yf-y2+1)+xinp10*(y2-yf);
%                 yinp1 = xinp01*(yf-y2+1)+xinp11*(y2-yf);
%                 out = yinp1*(z2-zf)+yinp0*(zf-z2+1);
%                 projection(sz,sy)=projection(sz,sy)+out;
%                 x=x+stepx;
%                 y=y+stepy;
%                 z=z+stepz;
%             end
%             startx=startx+steppixx;
%             starty=starty+steppixy;
%         end
%         startz=startz+steppixz;
%         if startz>105
%             break;
%         end
%     end
end
        