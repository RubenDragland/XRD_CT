%%SAVE_MASK_NEXUS saves a mask in an h5 file
% To be compatible with the rest of the code, make sure that the
% mask is saved without any modifications to its orientation, i.e. data and 
% mask should be loaded with orientation [0 0], as shown in the example, 
% processed and then saved.
% 
% ** filename       filename of the HDF5 file
% ** mask           2D/3D array
%
% EXAMPLES:
%       %% load already existing mask, modify it and save it again
%
%       mask = io.nexus_read('./oldMask.h5', 'filter', 'mask');
%       mask(250:350,2) = 0;
%       io.save_mask_nexus('./newMask.h5', mask);
%
%
%       %% load data and mask; modify existing mask with create_mask
%       %% and save the new mask to disk
%
%       [data, mask] = io.nexus_read('~/Data10', 'scanNr', 250, 'filter', 'pilatus_1', 'mask', './oldMask.h5', 'orientation', [0 0]);
%       plotting.imagesc3D(data); axis equal tight xy;
%       beamline.create_mask('mask', mask); % save_mask_nexus is called from beamline.create_mask
%
%
%       %% load data and create new mask with create_mask
%
%       data = io.nexus_read('~/Data10', 'scanNr', 250, 'filter', 'pilatus_1', 'orientation', [0 0]);
%       plotting.imagesc3D(data); axis equal tight xy;
%       beamline.create_mask();       
%
%
% see also: io.nexus_read


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


function save_mask_nexus(filename, mask)


s = [];
s.entry.collection.mask = mask;

% make sure that we don't have relative paths
if strcmpi(filename(1), '~')
    filename = utils.abspath(filename);
end

io.HDF.save2hdf5(filename, s, 'overwrite', true, 'comp', 2);


end