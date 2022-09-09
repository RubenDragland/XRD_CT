addpath ..

% show_sSAXS_pixel_scattering currently uses the output of
% analyze_one_scan.m, typically called from SAXS_caller_template....m

% directory and file-name with fourier.mat file, you must set at least this
% field for using this template
dir_file_fourier = '/das/work/p11/p11548/analysis/fourier_components/04_SAXS/C2_AAA_De_and_BGR_and_air/03357_fourier.mat';

show_sSAXS_pixel_scattering(dir_file_fourier, ...
    'which_qrange_to_plot', 1, ... % number of the entry in the list of pixel/q ranges defined in SAXS_caller_template
    'fig_type','fig_asym_sym_int', ... % name of the plot to use, see show_sSAXS_pixel_scattering or SAXS_caller_template for the fmt field-names that can be specified here
    'show_diffraction_pattern',0, ... % 1 to display the raw data as well
    'compile_x12sa_filename_args', {}, ... % Arguments to the function utils.compile_x12sa_filename to find the name of the scatterin pattern frame
    'image_show_args', {}, ... % Arguments to the function plotting.image_show 
    'basedir_integrated_data', [], ... % If empty it is taken from fnames.basedir_integrated_data
    'plot_radial_integ_args', ... % parameters for plotting.plot_radial_integ
    {'QMulPow',0.0, ... % especially in low q it may improve visibility of peaks if I*q^1 up to I*q^4 is plotted.
    'PlotQ',1 ... % 1 to plot as a functgion of q rather than pixels
    });

% The above settings are good for the beamline. If using in your own server
% or in Ra you will need to add a few more paths and info to the function
% in order for it to find your data. See an example below, note you will
% still need to update with your own paths and your own eaccount and
% pgroup.

% show_sSAXS_pixel_scattering(dir_file_fourier, ...
%     'which_qrange_to_plot', 2, ... % number of the entry in the list of pixel/q ranges defined in SAXS_caller_template
%     'fig_type','fig_asym_sym_int', ... % name of the plot to use, see show_sSAXS_pixel_scattering or SAXS_caller_template for the fmt field-names that can be specified here
%     'show_diffraction_pattern',1, ... % 1 to display the raw data as well
%     'compile_x12sa_filename_args', {'BasePath','/das/work/p17/p17055/pilatus_1/', 'BaseName','e17055_'}, ... % Arguments to the function utils.compile_x12sa_filename to find the name of the scatterin pattern frame
%     'image_show_args', {}, ... % Arguments to the function plotting.image_show 
%     'basedir_integrated_data', '/das/work/p17/p17055/analysis/radial_integration/', ... % If empty it is taken from fnames.basedir_integrated_data
%     'plot_radial_integ_args', ... % parameters for plotting.plot_radial_integ
%     {'QMulPow',0.0, ... % especially in low q it may improve visibility of peaks if I*q^1 up to I*q^4 is plotted.
%     'PlotQ',1 ... % 1 to plot as a functgion of q rather than pixels
%     });


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

