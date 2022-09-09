%CONVERT_LOCATION
% Translate a given pixel location from one detector orientation to another
% 
% ** pos            pixel location (Y/X)
% ** framesz        size of the data frame
% ** vecIn          orientation vector of the input data frame
% ** vecTypeIn      orientation vector type of the input data frame;
%                   default: 'transpose_fliplr_flipud'
% ** vecout         orientation vector of the output data frame; default [0 0 0]
% ** vecTypeOut     orientation vector type of the output data frame;
%                   default: 'transpose_fliplr_flipud'
%
% EXAMPLE:
%       math.convert_location('pos', [450 520], 'framesz', [1000 1000], 'vecIn', [1 1 1])
%
% see also: math.applyTransform.m


%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2019 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%*-----------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing language this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS software package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.

function posOut = convert_location(varargin)

par = inputParser;
par.addParameter('pos', [], @isnumeric)
par.addParameter('framesz', [], @isnumeric)
par.addParameter('vecIn', [], @isnumeric)
par.addParameter('vecTypeIn', 'transpose_fliplr_flipud', @ischar)
par.addParameter('vecOut', [0 0], @isnumeric)
par.addParameter('vecTypeOut', 'transpose_fliplr_flipud', @ischar)

par.parse(varargin{:})
vars = par.Results;

% input checks
if length(vars.pos)~=2
    error('Y/X position required. Please check your input argument "pos".')
end
if length(vars.framesz)<2
    error('Frame size should be 2D.')
end

% convert input transform vector to transpose_fliplr_flipud
if ~strcmpi(vars.vecTypeIn, 'transpose_fliplr_flipud')
    vars.vecIn = math.convert_transform(vars.vecIn, vars.vecTypeIn, 'transpose_fliplr_flipud');
end


% input
posOut = vars.pos;
if vars.vecIn(3)
    posOut(1) = vars.framesz(1) - vars.pos(1) + 1;
end
if vars.vecIn(2)
    posOut(2) = vars.framesz(2) - vars.pos(2) + 1;
end
if vars.vecIn(1)
    posOut = fliplr(posOut);
end

if iscolumn(posOut)
    posOut = posOut';
end


% output

if any(vars.vecOut~=0)
    % get new frame size
    if vars.vecIn(1)
        vars.framesz = fliplr(vars.framesz);
    end
    vars.pos = posOut;
    % convert output transform vector to transpose_fliplr_flipud
    if ~strcmpi(vars.vecTypeOut, 'transpose_fliplr_flipud')
        vars.vecOut = math.convert_transform(vars.vecOut, vars.vecTypeOut, 'transpose_fliplr_flipud');
    end
    if vars.vecOut(3)
        posOut(1) = vars.framesz(1) - vars.pos(1) + 1;
    end
    if vars.vecOut(2)
        posOut(2) = vars.framesz(2) - vars.pos(2) + 1;
    end
    if vars.vecOut(1)
        posOut = fliplr(posOut);
    end 


    
    if iscolumn(posOut)
        posOut = posOut';
    end

end
