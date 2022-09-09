%fluo_raw_to_spectra    Receives raw data structure from Falcon detector and returns array with arranged hyperspectra
%
% ** fluo_raw_Data         file path or path to a parent directory
%
% ** dead_time_correction   option (true or false) for not correcting for the dead time
%                           the default is to correct (true)
%
% returns:
% ++ fluo_data              requested fluorescence data
%
% EXAMPLES:
%       % specify the full path and nexus as data type
%       falcon = io.nexus_read('~/Data10','scanNr',scans(ii),'filter','falcon');
%       falcon_spectra = squeeze(fluo.fluo_raw_to_spectra(falcon));  Example for one detector

%
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

function fluo_data = fluo_raw_to_spectra(fluo_raw_data_str, varargin)

par = inputParser;
par.addParameter('dead_time_correction', true, @islogical)
par.parse(varargin{:})
vars = par.Results;

dead_time_correction = vars.dead_time_correction;

number_of_spectra = fluo_raw_data_str.CurrentPixel(end);
fluo_data = fluo_read_loc(fluo_raw_data_str, dead_time_correction, number_of_spectra);

return
end




function fluo_data_loc = fluo_read_loc(fluo_data_loc, dead_time_correction, number_of_spectra)
if isstruct(fluo_data_loc)
    try
        try
            buffer_data = fluo_data_loc.entry.data.data;
        catch ME
        end
        try
            buffer_data = fluo_data_loc.data;
        catch
        end
    catch ME
        error('Not implemented fluorescence structure. Falcon decode buffer accepts /entry/data/data or already extracted array.')
    end
else
    buffer_data = fluo_data_loc;
end

decoded_data = fluo.falcon_decode_buffers(buffer_data);
fluo_data_loc= decoded_data.Data(:, :, 1:number_of_spectra);

if dead_time_correction == 1
    live_time = decoded_data.LiveTime(1:number_of_spectra);
    real_time = decoded_data.RealTime(1:number_of_spectra);
    correction_coef = real_time ./ live_time;
    
    for ch_num = 1:number_of_spectra
        fluo_data_loc(:, ch_num) = fluo_data_loc(:, ch_num) * correction_coef(ch_num);
    end
end
end
