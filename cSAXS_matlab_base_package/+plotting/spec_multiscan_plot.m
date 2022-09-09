% function h = spec_multiscan_plot(spec_dat_file, scans, ['counter',counter_name])
%
%SPEC_MESH_PLOT Plot the counter values from a mesh scan. In the future to
%   be generalized to other types of scans
%
% inputs
%
% ** spec_dat_file  Filename and path of the SPEC dat file. A base path can
%                   also be given and the code will search for the file
% ** scans          Scan numbers
%
%*optional*
% ** 'counter',counter_name     counter_name is the name of the SPEC
%                               counter, by default it will be 'bpm4i'
%
% returns
% ++ h              Figure handle
%
%
% see also: 
%
% EXAMPLES:
%   plotting.spec_multiscan_plot('~/Data10',[263:268])
%   plotting.spec_multiscan_plot('~/Data10',[263:268], 'counter','diode')

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


function h = spec_multiscan_plot(spec_dat_file, scans, varargin)

% spec_dat_file = '~/Data10';
% scans = [263:268];


% set default values for the variable input arguments:
counter = 'bpm4i';

% check minimum number of input arguments
if (nargin < 2)
    error('At least the spec dat filename and scan number have to be specified as input parameters.');
end

% accept cell array with name/value pairs as well
no_of_in_arg = nargin;

% check number of input arguments
if (rem(no_of_in_arg,2) ~= 0)
    error('The optional parameters have to be specified as ''name'',''value'' pairs');
end


% parse the variable input arguments
vararg_remain = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch lower(name)
        case 'counter' 
            counter = value;
        otherwise
            vararg_remain{end+1} = name; %#ok<AGROW>
            vararg_remain{end+1} = value; %#ok<AGROW>
    end
end

Sp = io.spec_read(spec_dat_file,'ScanNr',scans);
if numel(Sp) == 1
    S{1} = Sp;
else
    S = Sp;
end
scanstring = strsplit(S{1}.S);
scantype = scanstring{3};

switch lower(scantype)
    case 'ascan'
        name_axis_fast = scanstring{4};
        N_fast = str2num(scanstring{7})+1;
    otherwise
        error(['Scan type ' lower(scantype) ' is not defined for this function'])
end
    

if nargout > 0
    h = figure(1);
else
    figure(1)
end

clf
hold on
legendtext = {};
for ii = 1:numel(scans)
    plot(S{ii}.(name_axis_fast),S{ii}.(counter))
    legendtext{ii} = sprintf('Scan %05d',scans(ii));
end
hold off
legend(legendtext)
xlabel(name_axis_fast)
ylabel(counter)

end








% data = io.spec_read('~/Data10/','ScanNr',scans);
% 
% figure(1);
% clf
% hold on
% for ii = 1:numel(scans)
%     plot(data{ii}.idgap,data{ii}.bpm4i)
% end
% hold off
% legend