% [ str ] = q_str_par( pixels, par, fmt )

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%*-----------------------------------------------------------------------*
%
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing la this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS scanning SAXS package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% Additionally, any publication using the package, or any translation of the 
%     code into another computing language should cite:
%    O. Bunk, M. Bech, T. H. Jensen, R. Feidenhans'l, T. Binderup, A. Menzel 
%    and F Pfeiffer, “Multimodal x-ray scatter imaging,” New J. Phys. 11,
%    123016 (2009). (doi: 10.1088/1367-2630/11/12/123016)
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ str ] = q_str_par( pixels, par, fmt )

if (nargin < 5)
    digits = 0;
end

ind = 1;
ind_i = ind;
if (isfield(fmt,'Angstrom'))
    nm_conversion = 1.0;
else
    nm_conversion = 0.1;
end
str = sprintf( fmt.q_format,4*pi / (par.lambda_A*nm_conversion) * sin(atan(pixels(ind)*par.pix_size/par.det_dist) /2) );
plot_ind = 1;
while (ind <= length(pixels))
    if ((ind == length(pixels))  || ...
        ((pixels(ind+1) ~= pixels(ind)+1) && (pixels(ind+1) ~= pixels(ind)-1)))
        if (ind ~= ind_i)
            str = [ str ' - ' ...
                sprintf(fmt.q_format,4*pi / (par.lambda_A*nm_conversion) * sin(atan(pixels(ind)*par.pix_size/par.det_dist) /2)) ];
            plot_ind = plot_ind +1;
        end
        if (ind < length(pixels))
            str = [ str ',  ' ];
            if (rem(plot_ind,5) == 0)
                str = [ str '\newline' ];
            end
            str = [ str sprintf(fmt.q_format,4*pi / (par.lambda_A*nm_conversion) * sin(atan(pixels(ind+1)*par.pix_size/par.det_dist) /2)) ];
            plot_ind = plot_ind +1;
            ind_i = ind+1;
        end
    end
    ind = ind +1;
end


if (isfield(fmt,'Angstrom'))
    str = [ str ' A^{-1}' ];
else
    str = [ str ' nm^{-1}' ];
end
