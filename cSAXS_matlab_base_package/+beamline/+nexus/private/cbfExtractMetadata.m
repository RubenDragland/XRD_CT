%CBFEXTRACTMETADATA read header and convert it into dataset/attribute
%structure for a nexus file
% 
% ** h              cbf header
%
% returns:
% ++ m              dataset/attribute structure
%
% see also: beamline.nexus.convert2NEXUS

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

function m = cbfExtractMetadata(h)
import beamline.nexus.*

[m.count_time.Value, m.count_time.Attributes.units] = getVal(h, '# Exposure_time', true);
[m.count_period.Value, m.count_period.Attributes.units] = getVal(h, '# Exposure_period', true);
[m.count_cutoff.Value, m.count_cutoff.Attributes.units] = getVal(h, '# Count_cutoff', true);
[m.tau.Value, m.tau.Attributes.units] = getVal(h, '# Tau', true);
[m.threshold_setting.Value, m.tau.Attributes.units] = getVal(h, '# Threshold_setting', true);
[m.gain, ~] = getVal(h, '# Gain_setting:', false);
[m.trim_file, ~] = getVal(h, '# Trim_file:', false);

end

function [v, unit] = getVal(h, tag, getUnit)
unit = [];
for ii=1:length(h)
    if contains(h{ii}, tag)
        res = strsplit(h{ii}, ' ');
        if getUnit
            jj=1;
            while jj<=length(res)
                tmp = str2double(res{jj});
                if isnumeric(tmp) && ~isnan(tmp)
                    v = tmp;
                    unit = res{jj+1};
                    break;
                end
                jj = jj+1;
            end
        else
            v = strjoin(res(3:end), ' ');
        end
    end
end
end
