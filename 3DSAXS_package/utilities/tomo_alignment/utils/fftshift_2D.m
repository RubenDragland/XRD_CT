function [x,idx]  = fftshift_2D(x, use_builtin)
% -----------------------------------------------------------------------
% This file is part of the PTYCHOMAT Toolbox
% Author: Michal Odstrcil, 2016
% License: Open Source under GPLv3
% Contact: ptychomat@gmail.com
% Website: https://bitbucket.org/michalodstrcil/ptychomat
% -----------------------------------------------------------------------
% Description:   faster version of matlab fftshift to work for stack of
% 2D images 

if nargin == 1 || ~use_builtin

    numDims = 2;
    idx = cell(1, numDims);
    for k = 1:numDims
        m = size(x, k);
        p = ceil(m/2);
        idx{k} = [p+1:m 1:p];
    end
    x = x(idx{:},:);
else
    x = builtin('fftshift', x);
end

end
