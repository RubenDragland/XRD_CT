%REDUCE_XRAY_DATABASE reduces the X-ray database to the elements that
%                     current beam energy can excite, either scan number or
%                     beam energy must be provided
%
% ** fls             fluo structure
%
% *optional*
% ** scanNr                 scan number
% ** beam_energy            beam energy in keV
%
% returns:
% ++ fls             fluo structure with reduced X-ray database
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

function [fls] = reduce_xray_database(fls, varargin)

par = inputParser;
par.addParameter('scanNr', -1, @isnumeric)
par.addParameter('beam_energy', -1, @isnumeric)

par.parse(varargin{:})
vars = par.Results;

if (vars.scanNr < 0) && (vars.beam_energy < 0)
    error('Either scan number or the beam energy are required.')
elseif (vars.scanNr > 0) && (vars.beam_energy > 0)
    warning(['No need to provide both the scan number and the beam energy\n' ...
             'Will use the manually entered beam_energy for the reduction of X-ray lines and edges database.'])
    vars.ScanNr = -1;
end

if vars.scanNr > 0
    if strcmp(fls.proc.file_type, 'nexus')
        beam_energy = io.nexus_read(fls.path.base_path, 'scanNr', vars.scanNr(1), 'filter', 'mokev');
    elseif strcmp(fls.proc.file_type, 'raw')
        motor_pos = io.spec_read(fls.path.base_path, 'ScanNr', vars.scanNr(1));
        beam_energy = motor_pos.mokev;
    end
elseif vars.beam_energy > 0
    beam_energy = vars.beam_energy;
else
    error('Either the scan number of the beam energy are wrong.')
end

indices = [fls.xray_database.edges.energy] - beam_energy <= 0;
edges = fls.xray_database.edges(indices);

indices = [fls.xray_database.lines.energy] - beam_energy <= 0;
lines = fls.xray_database.lines(indices);

reduced_lines = [];
prev_element = edges(1).element;
indices = [lines.element] == edges(1).element;
temp_lines = lines(indices);
if ~isempty(temp_lines)
    sub_indices = [temp_lines.initial_level] == edges(1).edge;
    reduced_lines = [reduced_lines; temp_lines(sub_indices)];
end

for i=2:length(edges)
    if edges(i).element == prev_element
        if ~isempty(temp_lines)
            sub_indices = [temp_lines.initial_level] == edges(i).edge;
            reduced_lines = [reduced_lines; temp_lines(sub_indices)];
        end
    elseif edges(i).element ~= prev_element
        prev_element = string(edges(i).element);
        indices = [lines.element] == edges(i).element;
        temp_lines = lines(indices);
        if ~isempty(temp_lines)
            sub_indices = [temp_lines.initial_level] == edges(i).edge;
            reduced_lines = [reduced_lines; temp_lines(sub_indices)];
        end
    end
end

fls.xray_database.reduced_edges = edges;
fls.xray_database.reduced_lines = reduced_lines;
end