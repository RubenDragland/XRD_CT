% analyze_one_scan(fnames,par,fmt,processing)
% Call function without arguments for a detailed explanation of its use

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


function [] = analyze_one_scan(fnames,par,fmt,processing)
import beamline.radial_integ

set(0, 'DefaultAxesfontsize', 10);
set(0, 'DefaultAxeslinewidth', 0.5, 'DefaultAxesfontsize', 10);
set(0, 'DefaultLinelinewidth', 0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setting default values for position of figures %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% position and size of figures
% fmt.fig_pos_sym_int          = [   6   579   569   366 ];
% fmt.fig_pos_deg_orient       = fmt.fig_pos_sym_int + [ 40 -40 0 0 ];
% fmt.fig_pos_orient           = fmt.fig_pos_deg_orient + [ 40 -40 0 0 ];
% fmt.fig_pos_comb_all         = [ 689   518   588   427 ];
% fmt.fig_pos_comb_asym        = fmt.fig_pos_comb_all + [ -40 -40 0 0 ];
% fmt.fig_pos_scatt_asym       = [ 882    81   374   269 ];
% fmt.fig_pos_cos_dev          = fmt.fig_pos_scatt_asym + [ -40 40 0 0 ];
% fmt.fig_pos_int_pow_law_exp  = fmt.fig_pos_comb_all + [ -40 40 0 0 ];
% fmt.fig_pos_fit_peak.pos     = fmt.fig_pos_comb_all + [ -80 40 0 0 ];
% fmt.fig_pos_fit_peak.ampl    = fmt.fig_pos_fit_peak.pos + [ -80 40 0 0 ];
% fmt.fig_pos_fit_peak.width   = fmt.fig_pos_fit_peak.ampl + [ -80 40 0 0 ];
% fmt.fig_pos_fit_peak.all     = fmt.fig_pos_fit_peak.width + [ -80 40 0 0 ];
% fmt.fig_pos_simple_peak = fmt.fig_pos_fit_peak;
%fig_pos_sym_asym         = [   6   126   468   819 ];
%fig_pos_orient           = [];
%fig_pos_comb_asym        = [ 855   467   421   478 ];
%fig_pos_comb_all         = [ 855   467   421   478 ];
%fig_pos_int_asym         = [ 902    61   374   669 ];
%fig_pos_int_pow_law_exp  = [ 515   489   450   456 ];

% display d-spacing in title and filename in this format
fmt.d_spacing_format = '%.3f';
fmt.d_spacing_file_format = '%05.2f';

% fmt.axis_on = 0;
% fmt.axis_xy = 0;
% fmt.fliplr = 0;
% fmt.flipud = 0;
% fmt.permutedim = 0;
fmt.cw_xoffs = 0.02;
% fmt.cw_yoffs = -0.01;
% fmt.title_on = 0;
% fmt.scale_bar_pos_pt = [ 20 145 ];
% fmt.scale_bar_x_mm = 1.00;
% fmt.scale_bar_y_mm = 0.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% directory for loading integrated data
if (isempty(fnames.basedir_integrated_data))
    % set default
    fnames.basedir_integrated_data = fullfile(fnames.data_dir, 'radial_integration/', fnames.rel_dir);
end
% directory and first part of filename for loading background data
if ((~isfield(fnames,'basename_bgr_data')) || (isempty(fnames.basename_bgr_data)))
    fnames.basename_bgr_data = fnames.basename_integrated_data;
end

% save Fourier components to this file
if (isempty(fnames.basedir_fourier))
    % set default
    fnames.basedir_fourier = fullfile(fnames.save_dir, 'fourier_components/', fnames.sub_dir, fnames.rel_dir );
end
if (isempty(fnames.filename_fourier))
    % set default
    fnames.filename_fourier = fullfile(fnames.basedir_fourier, [num2str(fnames.first_scan_no,'%05d') '_fourier.mat']);
end
% in case its not specified, it wil consider that snake_scan = 0
if ~isfield(par,'snake_scan')
    par.snake_scan = 0;
end
if ~isfield(par,'prepare_SASTT')
    par.prepare_SASTT = 0;
end
if ~isfield(par,'which_qrange')
    par.which_qrange = [];
end
if ~isfield(par,'fast_axis_x')
    par.fast_axis_x = 0;
end

if (par.prepare_SASTT)&&(~processing.load_data)&&(~isempty(par.qresolved_q))
    error(['You indicated you wanted to prepare SASTT data (par.prepare_SASTT) with q-resolution (par.qresolved_q), for these options you must reload the radial integrated data. Set processing.load_data = 1 and have fun.'])
end

% Check par.rsum
for ii = 1:numel(par.r_sum)
    if isempty(par.r_sum{ii})
        error(sprintf('par.r_sum{%d} is empty',ii));
    end
end

% ensure that the output directory for Fourier components exists
[mkdir_stat,mkdir_message] = mkdir(fnames.basedir_fourier);
if (~mkdir_stat)
    error('invalid directory %s: %s',fnames.basedir_fourier,mkdir_message);
end
if ((mkdir_stat) && (isempty(mkdir_message)))
    fprintf('The output directory %s for the Fourier components has been created.\n',...
        fnames.basedir_fourier);
end
if (processing.save_fig)
    % save figures to this directory
    if (isempty(fnames.basedir_figures))
        % set default
        fnames.basedir_figures = ...
            [ fnames.save_dir 'figures_saxs/' ];
    end
    fnames.basedir_figures = [ fnames.basedir_figures fnames.sub_dir fnames.rel_dir ];
    
    % ensure that the output directory for figures exists
    [mkdir_stat,mkdir_message] = mkdir(fnames.basedir_figures);
    if (~mkdir_stat)
        error('invalid directory %s: %s',fnames.basedir_figures,mkdir_message);
    end
    if ((mkdir_stat) && (isempty(mkdir_message)))
        fprintf('The output directory %s for the figures has been created.\n',...
            fnames.basedir_figures);
    end
    
    subdir_name = [ fnames.basedir_figures '/eps/'];
    [mkdir_stat,mkdir_message] = mkdir(subdir_name);
    if ((mkdir_stat) && (isempty(mkdir_message)))
        fprintf('The output directory %s for the figures has been created.\n',...
            subdir_name);
    end
    
    subdir_name = [ fnames.basedir_figures '/jpg/'];
    [mkdir_stat,mkdir_message] = mkdir(subdir_name);
    if ((mkdir_stat) && (isempty(mkdir_message)))
        fprintf('The output directory %s for the figures has been created.\n',...
            subdir_name);
    end
    
    subdir_name = [ fnames.basedir_figures '/fig/'];
    [mkdir_stat,mkdir_message] = mkdir(subdir_name);
    if ((mkdir_stat) && (isempty(mkdir_message)))
        fprintf('The output directory %s for the figures has been created.\n',...
            subdir_name);
    end
end

% Creating the title string
fmt.title_base = strrep(fnames.rel_dir,'/','');
fmt.title_base = [ strrep(fmt.title_base,'_','\_') ': ' ];

if ((~isfield(fmt,'title_user_start')) || (isempty(fmt.title_user_start)))
    fmt.title_user_start = '';
end

% load data from previously performed analysis
if (~processing.load_data)
    fid = fopen(fnames.filename_fourier);
    if (fid >= 0)
        fclose(fid);
        
        four_dat = [];
        
        % The load-data flag is used to force re-loading and re-processing.
        % It is not useful to load this from the previous processing run,
        % but the processing structure may be stored anyway.
        fprintf('loading %s\n',fnames.filename_fourier);
        load(fnames.filename_fourier);
        
        % put q-vector and detector distance in place
        if (isfield(four_dat.par,'q'))
            par.q = four_dat.par.q;
        end
        if (isfield(four_dat.par,'phi_det'))
            par.phi_det = four_dat.par.phi_det;
        end
        
        % enforce reprocessing if the integration ranges changed
        if (isfield(four_dat,'par'))
            if ((length(par.r_sum) ~= length(four_dat.par.r_sum)) || ...
                    (length(par.r_bgr) ~= length(four_dat.par.r_bgr)) || ...
                    (length(par.I_bgr_pow) ~= length(four_dat.par.I_bgr_pow)) || ...
                    (par.ind_max_x ~= four_dat.par.ind_max_x) || ...
                    (par.ind_max_y ~= four_dat.par.ind_max_y))
                processing.load_data = 1;
            end
            if (~processing.load_data)
                if ((nnz(par.r_bgr ~= four_dat.par.r_bgr)) || ...
                        (nnz(par.I_bgr_pow ~= four_dat.par.I_bgr_pow)))
                    processing.load_data = 1;
                end
            end
            if (~processing.load_data)
                for cmp_ind = 1:length(four_dat.par.r_sum)
                    if (((nnz(size(par.r_sum{cmp_ind}) ~= size(four_dat.par.r_sum{cmp_ind}))) || ...
                            (nnz(par.r_sum{cmp_ind} ~= four_dat.par.r_sum{cmp_ind}))))
                        processing.load_data = 1;
                        break;
                    end
                end
            end
            % check if peak-fitting parameters changed
            if (~processing.load_data)
                if ((any(size(par.r_sum) ~= size(four_dat.par.r_sum))) || ...
                    (any(par.peakfit ~= four_dat.par.peakfit)))
                    processing.load_data = 1;
                end
            end
        end
        if (~processing.load_data)
            par.ind_max_y = size(four_dat.f1_amp,1) -1;
            par.ind_max_x = size(four_dat.f1_amp,2) -1;
        else
            fprintf('Mismatch found, reprocessing ...\n');
        end
    end
end



if ((~exist('four_dat','var')) || (processing.load_data))
    reload_data = 1;
else
    reload_data = 0;
end

% perform Fourier analysis
if (reload_data)
    [four_dat, data_process] = fourier_analysis_linewise(fnames,par,fmt);
    fprintf('saving data to %s\n',fnames.filename_fourier);
    save(fnames.filename_fourier, 'four_dat',  'data_process', '-v6');
end
% save parameters in return structure to avoid that current settings are
% overwritten upon reloading the result of previous processing
four_dat.par = par;
four_dat.fmt = fmt;
four_dat.fnames = fnames;
plot_fourier_result(four_dat, processing);

if (par.prepare_SASTT)
    prepare_SASTT(four_dat, data_process);
end
