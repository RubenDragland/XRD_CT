%function plot_fourier_result(four_dat, processing)

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

function plot_fourier_result(four_dat, processing)

fnames = four_dat.fnames;
par = four_dat.par;
fmt = four_dat.fmt;

if isfield(four_dat,'integ')
    integ = four_dat.integ;
else
    warning('four_dat.integ does not exist, using legacy compatibility, checking in four_dat.par for phi_det')
    integ.phi_det = four_dat.par.phi_det;
    integ.q       = four_dat.par.q;
end

if ~isfield(processing,'which_qrange_to_plot')
    processing.which_qrange_to_plot = [];
end

fmt.printer_name = 'wsla_x12sa';
fmt.x = (0:par.ind_max_x) * par.x_scale;
fmt.y = (0:par.ind_max_y) * par.y_scale;

fmt.axes_font_size = 10;
fmt.annotation_font_size = 10;
fmt.colorbar_font_size = 10;

% start position and position offset for each figure
xs =  5;
ys = 50;
xo = 55;
yo = 30;

% set default figure size and position
if (~isfield(fmt,'fig_pos_sym_int')) || (isempty(fmt.fig_pos_sym_int))
    fmt.fig_pos_sym_int =  [ xs+01*xo ys+01*yo 700 600 ];
end
if (~isfield(fmt,'fig_pos_deg_orient')) || (isempty(fmt.fig_pos_deg_orient))
    fmt.fig_pos_deg_orient =  [ xs+02*xo ys+02*yo 700 600 ];
end
if (~isfield(fmt,'fig_pos_orient')) || (isempty(fmt.fig_pos_orient))
    fmt.fig_pos_orient =  [ xs+03*xo ys+03*yo 700 600 ];
end
if (~isfield(fmt,'fig_pos_comb_asym')) || (isempty(fmt.fig_pos_comb_asym))
    fmt.fig_pos_comb_asym =  [ xs+04*xo ys+04*yo 700 600 ];
end
if (~isfield(fmt,'fig_pos_comb_all')) || (isempty(fmt.fig_pos_comb_all))
    fmt.fig_pos_comb_all =  [ xs+05*xo ys+05*yo 700 600 ];
end
if (~isfield(fmt,'fig_pos_scatt_asym')) || (isempty(fmt.fig_pos_scatt_asym))
    fmt.fig_pos_scatt_asym =  [ xs+20*xo ys+20*yo 700 600 ];
end
if (~isfield(fmt,'fig_pos_cos_dev')) || (isempty(fmt.fig_pos_cos_dev))
    fmt.fig_pos_cos_dev =  [ xs+21*xo ys+21*yo 700 600 ];
end
if (~isfield(fmt,'fig_pos_pow_law_exp')) || (isempty(fmt.fig_pos_int_pow_law_exp))
    fmt.fig_pos_int_pow_law_exp = [ xs+06*xo ys+06*yo 700 600 ];
end

if (~isfield(fmt,'fig_pos_fit_peak'))
    fmt.fig_pos_fit_peak = {};
end
if (~isfield(fmt,'fig_pos_fit_peak.pos') || (isempty(fmt.fig_pos_fit_peak.pos)))
    fmt.fig_pos_fit_peak.pos = [ xs+07*xo ys+07*yo 700 600 ];
end
if (~isfield(fmt,'fig_pos_fit_peak.ampl') || (isempty(fmt.fig_pos_fit_peak.ampl)))
    fmt.fig_pos_fit_peak.ampl = [ xs+08*xo ys+08*yo 700 600 ];
end
if (~isfield(fmt,'fig_pos_fit_peak.width') || (isempty(fmt.fig_pos_fit_peak.width)))
    fmt.fig_pos_fit_peak.width = [ xs+09*xo ys+09*yo 700 600 ];
end
if (~isfield(fmt,'fig_pos_fit_peak.all') || (isempty(fmt.fig_pos_fit_peak.all)))
    fmt.fig_pos_fit_peak.all = [ xs+10*xo ys+10*yo 700 600 ];
end

if (~isfield(fmt,'fig_pos_simple_peak'))
    fmt.fig_pos_simple_peak = {};
end
if (~isfield(fmt,'fig_pos_simple_peak.pos') || (isempty(fmt.fig_pos_simple_peak.pos)))
    fmt.fig_pos_simple_peak.pos = [ xs+11*xo ys+11*yo 700 600 ];
end
if (~isfield(fmt,'fig_pos_simple_peak.ampl') || (isempty(fmt.fig_pos_simple_peak.ampl)))
    fmt.fig_pos_simple_peak.ampl = [ xs+12*xo ys+12*yo 700 600 ];
end
if (~isfield(fmt,'fig_pos_simple_peak.width') || (isempty(fmt.fig_pos_simple_peak.width)))
    fmt.fig_pos_simple_peak.width = [ xs+13*xo ys+13*yo 700 600 ];
end
if (~isfield(fmt,'fig_pos_simple_peak.all') || (isempty(fmt.fig_pos_simple_peak.all)))
    fmt.fig_pos_simple_peak.all = [ xs+14*xo ys+14*yo 700 600 ];
end
if (~isfield(fmt,'parallel_plotting'))
    fmt.parallel_plotting = true;
end

% colow wheel radius
cw_rad = 50;

% standard gray-scale color map
fmt.cm_gray = gray(512);

% interpolation of the plotted data, 0 means off, greater than 1 will
% become slow
fmt.interp_to = fmt.interpol;

if (processing.movie_fig)
    filename_movie = fullfile(fnames.basedir_figures, [par.sample_name '_movie.avi']);
    fprintf('writing frames to %s\n',filename_movie);
    movie_obj1 = VideoWriter(filename_movie);
    movie_obj1.FrameRate = 1;
    open(movie_obj1);
    movie_figure = figure(98);
else
    movie_figure = [];
    movie_obj1 = [];
end


% calculate a color wheel
cw_size = 2*cw_rad +1;
[ cw_x, cw_y ] = meshgrid( (1:cw_size)-cw_rad-0.5, (1:cw_size)-cw_rad-0.5 );    
[ theta, rho ] = cart2pol( cw_x, cw_y );
theta = (-theta+pi/2) / pi;
% theta = theta  / pi;
ind_wrap = find(theta < 0);
theta(ind_wrap) = theta(ind_wrap) +1;
ind_wrap = find(theta >= 1);
theta(ind_wrap) = theta(ind_wrap) -1;
[ind1,ind2] = find(rho < cw_rad);
cw_ind_2d = sub2ind(size(theta),ind1,ind2);
% choose a white background for the colour wheel
cw_hsv = ones(cw_size,cw_size,3);
cw_hsv(:,:,2) = 0.0;
%
cw_ind1 = sub2ind(size(cw_hsv),ind1,ind2,  ones(size(ind1)));
cw_hsv(cw_ind1) = theta(cw_ind_2d);
cw_ind2 = sub2ind(size(cw_hsv),ind1,ind2,2*ones(size(ind1)));
cw_hsv(cw_ind2) = 0.999;
cw_ind3 = sub2ind(size(cw_hsv),ind1,ind2,3*ones(size(ind1)));
cw_hsv(cw_ind3) = rho(cw_ind_2d) / (cw_rad);
cw_rgb = hsv2rgb(cw_hsv);
% flip color wheel according to the data
% if ((isfield(fmt,'fliplr')) && (fmt.fliplr))
%     cw_rgb = flipdim(cw_rgb,2);
% end
% if ((isfield(fmt,'flipud')) && (fmt.flipud))
%     cw_rgb = flipdim(cw_rgb,1);
% end
if ((isfield(fmt,'permutedim')) && (fmt.permutedim))
    dim_order = [ 2 1 3 ];
    cw_rgb = permute(cw_rgb,dim_order);
end

            
ind_max = length(par.r_sum);


if isempty(processing.which_qrange_to_plot)
    indices = 1:ind_max;
else
    indices = processing.which_qrange_to_plot;
end

% create temporary varibles to make the slicing obvious to Matlab
% four_dat
f1_amp = four_dat.f1_amp;
f2_amp = four_dat.f2_amp;
f2_phase = four_dat.f2_phase;
degree_orient = four_dat.degree_orient;
orientation = four_dat.orientation;
if (isfield(four_dat,'I_dev'))
    I_dev = four_dat.I_dev;
else
    I_dev = [];
end
if (isfield(four_dat,'I_cos_dev'))
    I_cos_dev = four_dat.I_cos_dev;
else
    I_cos_dev = [];
end
if (isfield(four_dat,'I_pow'))
    I_pow = four_dat.I_pow;
else
    I_pow = [];
end
if (isfield(four_dat,'fit_peak'))
    fit_peak = four_dat.fit_peak;
else
    fit_peak = [];
end
if (isfield(four_dat,'simple_peak'))
    simple_peak = four_dat.simple_peak;
else
    simple_peak = [];
end

% resulting axes handles
axes_handles = cell(max(indices),1);

% create the specified plots for each q-range
% use for instead of parfor to see the plots on the screen rather than just
% having them created in a buffer and stored in different file formats
if fmt.parallel_plotting
    parfor ind=indices
        precomputation_and_plotting(ind,integ,par,fnames,fmt,processing,f1_amp,f2_amp,f2_phase,degree_orient,I_dev,I_cos_dev,I_pow,orientation,fit_peak, simple_peak,cw_rgb)
    end
else
    for ind=indices
        precomputation_and_plotting(ind,integ,par,fnames,fmt,processing,f1_amp,f2_amp,f2_phase,degree_orient,I_dev,I_cos_dev,I_pow,orientation,fit_peak, simple_peak,cw_rgb)
    end    
end
   

if (processing.movie_fig)
    close(movie_obj1);
end

% determine number of axes handles and transfer to normal array and connect
% all plots to zooming together
no_of_axes_handles = 0;
for ind=indices
    no_of_axes_handles = no_of_axes_handles + length(axes_handles{ind});
end
axes_handle_vec = zeros(no_of_axes_handles,1);
if (no_of_axes_handles > 0)
    no_of_axes_handles = 0;
    for ind=indices
        if (~isempty(axes_handles{ind}))
            axes_handle_vec((no_of_axes_handles+1):(no_of_axes_handles + length(axes_handles{ind}))) = axes_handles{ind};
            no_of_axes_handles = no_of_axes_handles + length(axes_handles{ind});
        end
    end

    % if the plots still exist (for example not in parallel mode where they
    % don't appear at all), then all plots zoom together 
    if (all(isgraphics(axes_handle_vec))) && (processing.plot_data)
        linkaxes(axes_handle_vec,'xy');
    end
end

% Reassign the modified structures to four_dat just in case it will be
% needed to be returned in the future
four_dat.fnames = fnames;
four_dat.par    = par;
four_dat.fmt    = fmt;
four_dat.integ  = integ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Start of auxiliary functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function precomputation_and_plotting(ind,integ,par,fnames,fmt,processing,f1_amp,f2_amp,f2_phase,degree_orient,I_dev,I_cos_dev,I_pow,orientation,fit_peak, simple_peak,cw_rgb)

Idev    = [];
Icosdev = [];
Ipow    = [];
f1amp = f1_amp(:,:,ind);
f2amp = f2_amp(:,:,ind);
d = degree_orient(:,:,ind);
if (~isempty(I_dev))
    Idev = I_dev(:,:,ind);
end
if (~isempty(I_cos_dev))
    Icosdev = I_cos_dev(:,:,ind);
end
if (~isempty(I_pow))
    Ipow = I_pow(:,:,ind);
end

if (isnan(f2_phase(:,:,ind)))
    fprintf('Skipping ind=%d due to phase=nan\n',ind);
    return
end
% phase scaled to range [0 1]
f2phase_01 = f2_phase(:,:,ind) - integ.phi_det(1)*pi/90; % 90 instead of 180 because of symmetric intensity calculation

% peak fitting: check if non-zero peak-position data are available,
% if yes, provide scaled arrays for plotting
if (~isempty(fit_peak)) && (sum(sum(fit_peak.pos(:,:,ind))) ~= 0)
    [ fitpeak ] = prepare_peak_fit_plot_data(fit_peak, ind, integ.q, fmt.std_scale);
else
    fitpeak = {};
end

% simple peak analysis: check if non-zero peak-position data are available,
% if yes, provide scaled arrays for plotting
if (~isempty(simple_peak)) && (sum(sum(simple_peak.pos(:,:,ind))) ~= 0)
    [ simplepeak ] = prepare_peak_fit_plot_data(simple_peak, ind, integ.q, fmt.std_scale);
else
    simplepeak = {};
end



ignore_f2p_sign = 0;

if (ignore_f2p_sign)
    f2phase_01  = f2phase_01 / pi;
else
    f2phase_01  = f2phase_01 / (2*pi) + 0.5;
end
ind_wrap = find(f2phase_01 < 0.0);
f2phase_01(ind_wrap) = f2phase_01(ind_wrap) + 1.0;
ind_wrap = find(f2phase_01 >= 1.0);
f2phase_01(ind_wrap) = f2phase_01(ind_wrap) - 1.0;

% degree of orientation, should be a value between 0 and 1
d_lin_med = reshape(medfilt2(d),numel(d),1);
d_std = std(d_lin_med);
d_med = median(d_lin_med);

% histogram of the degree of orientation
d_hist_bins = 0.000:0.001:1.000;
d_hist_data = histc(d(:),d_hist_bins);
d_hist = { d_hist_bins; d_hist_data };


f1_lin_med = reshape(medfilt2(f1amp),size(f1amp,1)*size(f1amp,2),1);
f1_std = std(f1_lin_med);
f1_med = median(f1_lin_med);

f2_lin_med = reshape(medfilt2(f2amp),size(f2amp,1)*size(f2amp,2),1);
f2_std = std(f2_lin_med);
f2_med = median(f2_lin_med);

if (~isempty(Idev))
    Idev_lin_med = reshape(medfilt2(Idev),numel(Idev),1);
    Idev_std = std(Idev_lin_med);
    Idev_med = median(Idev_lin_med);
else
    Idev_med= [];
    Idev_std = [];
end

if (~isempty(Icosdev))
    Icosdev_lin_med = reshape(medfilt2(Icosdev),numel(Icosdev),1);
    Icosdev_std = std(Icosdev_lin_med);
    Icosdev_med = median(Icosdev_lin_med);
else
    Icosdev_std = [];
    Icosdev_med = [];
end

if (~isempty(Ipow))
    Ipow_lin_med = reshape(medfilt2(Ipow),numel(Ipow),1);
    Ipow_std = std(Ipow_lin_med);
    Ipow_med = median(Ipow_lin_med);
else
    Ipow_med = [];
    Ipow_std = [];
end

% combined plot all:
pl_comb_all = zeros(size(f2phase_01,1),size(f2phase_01,2),3);

% scale the hue with the orientation
pl_comb_all(:,:,1) = f2phase_01;

% scale the saturation with the asymmetric amplitude
if (isfield(fmt,'asymmetric_amplitude_scale_max')) && (~isempty(fmt.asymmetric_amplitude_scale_max))
    s_scale = fmt.asymmetric_amplitude_scale_max;
else
    s_scale = (f2_med + fmt.std_scale * f2_std);
end
f2_scaled = 1.0 / s_scale * f2amp;
f2_scaled( f2_scaled > 1.0 ) = 1.0;
pl_comb_all(:,:,2) = f2_scaled;


% scale the value with the total amplitude
famp = f1amp;
f_lin = reshape(medfilt2(famp),size(famp,1)*size(famp,2),1);
f_std = std(f_lin);
f_med = median(f_lin);
if (isfield(fmt,'symmetric_amplitude_scale_max')) && (~isempty(fmt.symmetric_amplitude_scale_max))
    v_scale = fmt.symmetric_amplitude_scale_max;
else
    v_scale = (f_med + fmt.std_scale * f_std);
end
f_scale = 1.0 / v_scale;
famp = f_scale * famp;
famp( famp > 1.0 ) = 1.0;
pl_comb_all(:,:,3) = famp;

pl_comb_all = hsv2rgb(pl_comb_all);


% combined plot asymmetric intensity:
pl_comb_asym = zeros(size(f2phase_01,1),size(f2phase_01,2),3);
% scale the hue with the orientation
pl_comb_asym(:,:,1) = f2phase_01;
% set the saturation to 1
pl_comb_asym(:,:,2) = 1.0;
% scale the value with the asymmetric amplitude
%     f_lin = reshape(medfilt2(f2amp),size(f2amp,1)*size(f2amp,2),1);
%     f_std = std(f_lin);
%     f_med = median(f_lin);
v_scale_asym = (f2_med + fmt.std_scale * f2_std);
f_scale = 1.0 / v_scale_asym;
f_scaled = f_scale * f2amp;
f_scaled( f_scaled > 1.0 ) = 1.0;
pl_comb_asym(:,:,3) = f_scaled;

pl_comb_asym = hsv2rgb(pl_comb_asym);


if ((processing.plot_data) || ...
        (processing.print_fig) || ...
        (processing.save_fig))
    
    ah = [];
    
    % symmetric intensity, i.e. amplitude of intensity component being
    % uniform across all azimuths
    if (fmt.fig_sym_int)
        ah= [ah plot_symmetric_amplitude(ind*10+101, f1amp,f1_med,f1_std, fnames,par,fmt,integ,processing,ind)];
    end
    
    % degree of orientation, i.e. asymmetric over total intensity
    % amplitude
    if (fmt.fig_orient_degree)
        ah= [ah plot_degree_of_orientation(ind*10+102, d,d_med,d_std, fnames,par,fmt,integ,processing,ind)];
    end
    
    % histogram of orientations (of the azimuthal angle with
    % highest intensity)
    if (fmt.fig_orient_degree_histogram)
        ah= [ah plot_degree_of_orientation_histogram(ind*10+103, d_hist,d_med,d_std, fnames,par,fmt,integ,processing,ind)];
    end
    
    % orientation of the intensity, i.e. azimuthal angle with
    % highest intensity
    if (fmt.fig_orientation)
        ah= [ah plot_orientation(ind*10+104, orientation(:,:,ind), ignore_f2p_sign,fnames,par,fmt,integ,processing,ind)];
    end
    
    % orientation and asymmetric intensity (i.e. intensity that
    % deviates from uniform azimuthal distribution)
    if (fmt.fig_asym_int)
        ah= [ah plot_orientation_and_asymmetric_intensity(ind*10+106,pl_comb_asym,d_med,d_std, fnames,par,fmt,integ,processing,cw_rgb,ind)];
    end
    
    % orientation symmetric and asymmetric intensity
    if (fmt.fig_asym_sym_int)
        ah= [ah plot_orientation_asymmetric_symmetric_intensity(ind*10+105, pl_comb_all,d_med,d_std, fnames,par,fmt,integ,processing,cw_rgb,ind)];
    end
    
    % exponent of intensity decay
    if ((fmt.fig_I_pow) && (exist('Ipow','var')))
        ah= [ah plot_intensity_decay_exponent(ind*10+109, Ipow,Ipow_med,Ipow_std, fnames,par,fmt,integ,processing,ind)];
    end
    
    % intensity deviation from point symmetry, i.e. from the first half
    % of the azimuthal segments exhibiting the same intensity as the
    % opposite segments in the second half.
    if ((fmt.fig_scatt_asym) && (exist('Idev','var')))
        ah= [ah plot_intensity_debviation(ind*10+107, Idev,Idev_med,Idev_std, fnames,par,fmt,integ,processing,ind)];
    end
    
    % deviation of the azimuthal intensity distribution from a cosine,
    % as the Fourier analysis works best for a sine/cosine like
    % distribution
    if ((fmt.fig_cos_dev) && (exist('Icosdev','var')))
        ah= [ah plot_absolute_intensity_deviation_from_cosine(ind*10+108, Icosdev,Icosdev_med,Icosdev_std, fnames,par,fmt,integ,processing,ind)];
    end
    
    % simple peak analysis plots
    if (~isempty(simplepeak)) && (isfield(fmt, 'fig_simple_peak'))
        ah((end+1):(end+4)) = plot_peak_fit_data(isfield(fmt.fig_simple_peak, 'pos') && (fmt.fig_simple_peak.pos), ...
            isfield(fmt.fig_simple_peak, 'ampl') && (fmt.fig_simple_peak.ampl), ...
            isfield(fmt.fig_simple_peak, 'width') && (fmt.fig_simple_peak.width), ...
            isfield(fmt.fig_simple_peak, 'all') && (fmt.fig_simple_peak.all), ...
            fmt.fig_pos_simple_peak, ...
            fnames,par,fmt,integ,processing, ...
            'simple-peak', ...
            simplepeak, ...
            ind);
    end
    
    % peak fitting plots
    if (~isempty(fitpeak)) && (isfield(fmt, 'fig_fit_peak'))
        ah((end+1):(end+4)) = plot_peak_fit_data(isfield(fmt.fig_fit_peak, 'pos') && (fmt.fig_fit_peak.pos), ...
            isfield(fmt.fig_fit_peak, 'ampl') && (fmt.fig_fit_peak.ampl), ...
            isfield(fmt.fig_fit_peak, 'width') && (fmt.fig_fit_peak.width), ...
            isfield(fmt.fig_fit_peak, 'all') && (fmt.fig_fit_peak.all), ...
            fmt.fig_pos_fit_peak, ...
            fnames,par,fmt,integ,processing, ...
            'fit-peak', ...
            fitpeak, ...
            ind);
    end
    
    axes_handles{ind} = ah;
end

% movie, i.e. for each q-range update a figure with sub-plots and store
% this as one frame of the movie file
if (processing.movie_fig)
    plot_movie_figure_and_store_frame(movie_figure,movie_obj1, ...
        f1amp,f1_med,f1_std, d,d_med,d_std, pl_comb_asym, pl_comb_all, ...
        fnames,par,fmt,integ,processing,ind);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ah] = plot_symmetric_amplitude(fig_no, f1amp,f1_med,f1_std, fnames,par,fmt,integ,processing,ind)
fmt.title_end = ', symmetric amplitude (cts/pixel)';
fmt.scale = [ f1_med-fmt.std_scale*f1_std f1_med+fmt.std_scale*f1_std ];
fmt.scale_min = 0;
fmt.scale_max = [];
fmt.colormap = fmt.cm_gray;
fmt.fig_pos = fmt.fig_pos_sym_int;
fnames.save_rel_name = 'symmetric_intensity';
ah = yet_another_plot(fig_no,fnames,par,fmt,integ,processing,f1amp,[],ind);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ah] = plot_degree_of_orientation(fig_no, degree_orient,d_med,d_std, fnames,par,fmt,integ,processing,ind)
fmt.title_end = ', degree of orientation';
fmt.scale = [ d_med-fmt.std_scale*d_std d_med+fmt.std_scale*d_std ];
fmt.scale_min = 0.0;
fmt.scale_max = 1.0;
fmt.colormap = fmt.cm_gray;
fmt.fig_pos = fmt.fig_pos_deg_orient;
fnames.save_rel_name = 'degree_of_orientation';
ah = yet_another_plot(fig_no,fnames,par,fmt,integ,processing,degree_orient,[],ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ah] = plot_degree_of_orientation_histogram(fig_no, d_hist,d_med,d_std, fnames,par,fmt,integ,processing,ind)
fmt.histplot = 1;
fmt.title_end = ', histogram of the degree of orientation';
fmt.scale = [ d_med-fmt.std_scale*d_std d_med+fmt.std_scale*d_std ];
fmt.scale_min = 0.0;
fmt.scale_max = 1.0;
fmt.colormap = fmt.cm_gray;
fmt.fig_pos = fmt.fig_pos_deg_orient;
fnames.save_rel_name = 'degree_of_orientation_histogram';
ah = yet_another_plot(fig_no,fnames,par,fmt,integ,processing,d_hist,[],ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ah] = plot_orientation(fig_no, orientation, ignore_f2p_sign,fnames,par,fmt,integ,processing,ind)
fmt.title_end = ', orientation (degrees)';
if (ignore_f2p_sign)
    fmt.scale = [0 90];
else
    fmt.scale = [-90 90];
end
fmt.scale_min = [];
fmt.scale_max = [];
fmt.colormap = flipud( hsv(361) );
fmt.fig_pos = fmt.fig_pos_orient;
fnames.save_rel_name = 'orientation';
ah = yet_another_plot(fig_no,fnames,par,fmt,integ,processing,orientation,[],ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ah] = plot_orientation_asymmetric_symmetric_intensity(fig_no, pl_comb_all,d_med,d_std, fnames,par,fmt,integ,processing,cw_rgb,ind)
fmt.title_end = ', orientation, asym. and sym. intensity';
fmt.scale = [ d_med-fmt.std_scale*d_std d_med+fmt.std_scale*d_std ];
fmt.scale_min = [];
fmt.scale_max = [];
fmt.colormap = [];
fmt.fig_pos = fmt.fig_pos_comb_all;
fnames.save_rel_name = 'orientation_asym_sym_int';
ah = yet_another_plot(fig_no,fnames,par,fmt,integ,processing,pl_comb_all,cw_rgb,ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ah] = plot_orientation_and_asymmetric_intensity(fig_no, pl_comb_asym,d_med,d_std, fnames,par,fmt,integ,processing,cw_rgb,ind)
fmt.title_end = ', orientation and asym. intensity';
fmt.scale = [ d_med-fmt.std_scale*d_std d_med+fmt.std_scale*d_std ];
fmt.scale_min = [];
fmt.scale_max = [];
fmt.colormap = [];
fmt.fig_pos = fmt.fig_pos_comb_asym;
fnames.save_rel_name = 'orientation_asym_int';
ah = yet_another_plot(fig_no,fnames,par,fmt,integ,processing,pl_comb_asym,cw_rgb,ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ah] = plot_intensity_debviation(fig_no, Idev,Idev_med,Idev_std, fnames,par,fmt,integ,processing,ind)
fmt.title_end = ', scattering asymmetry';
fmt.scale = [ Idev_med-fmt.std_scale*Idev_std Idev_med+fmt.std_scale*Idev_std ];
fmt.scale_min = [];
fmt.scale_max = [];
fmt.colormap = fmt.cm_gray;
fmt.fig_pos = fmt.fig_pos_scatt_asym;
fnames.save_rel_name = 'scattering_asymmetry';
ah = yet_another_plot(fig_no,fnames,par,fmt,integ,processing,Idev,[],ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ah] = plot_absolute_intensity_deviation_from_cosine(fig_no, Icosdev,Icosdev_med,Icosdev_std, fnames,par,fmt,integ,processing,ind)
fmt.title_end = ', abs. int. dev. from cos.';
fmt.scale = [ Icosdev_med-fmt.std_scale*Icosdev_std Icosdev_med+fmt.std_scale*Icosdev_std ];
fmt.scale_min = 0.0;
fmt.scale_max = [];
fmt.colormap = fmt.cm_gray;
fmt.fig_pos = fmt.fig_pos_cos_dev;
fnames.save_rel_name = 'abs_int_dev_from_cos';
ah = yet_another_plot(fig_no,fnames,par,fmt,integ,processing,Icosdev,ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ah] = plot_intensity_decay_exponent(fig_no, Ipow,Ipow_med,Ipow_std, fnames,par,fmt,integ,processing,ind)
fmt.title_end = ', intensity exponent';
fmt.scale = [ Ipow_med-fmt.std_scale*Ipow_std Ipow_med+fmt.std_scale*Ipow_std ];
fmt.scale_min = -5.0;
fmt.scale_max =  1.0;
fmt.colormap = fmt.cm_gray;
fmt.fig_pos = fmt.fig_pos_int_pow_law_exp;
fnames.save_rel_name = 'intensity_exponent';
ah = yet_another_plot(fig_no,fnames,par,fmt,integ,processing,Ipow,[],ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_movie_figure_and_store_frame(movie_figure,movie_obj1, ...
            f1amp,f1_med,f1_std, d,d_med,d_std, pl_comb_asym, pl_comb_all, ...
            fnames,par,fmt,integ,processing,ind)
% Currently movie_fig does not work with save_fig, print_fig.
% movie with radius display
figure(movie_figure);
clf
hold off;
set(gcf,'Name','movie frames 1');
set(gcf,'Position',[400   300   782   507]);
% background color
set(gcf,'Color','white');

fmt.subplot = 1;
fmt.title_on = 1;
fmt.d_spacing_center_pix = 1;
fig_no = ind*10+101;

subplot(2,2,1);       
fmt.title_end = ', symmetric amplitude';
fmt.scale = [ f1_med-fmt.std_scale*f1_std f1_med+fmt.std_scale*f1_std ];
fmt.scale_min = 0;
fmt.scale_max = [];
fmt.colormap = fmt.cm_gray;
fmt.fig_pos = fmt.fig_pos_sym_int;
fnames.save_rel_name = 'symmetric_intensity';
yet_another_plot(fig_no,fnames,par,fmt,integ,processing,f1amp,[],ind);


subplot(2,2,2);       
fmt.title_end = 'degree of orientation';
fmt.scale = [ d_med-fmt.std_scale*d_std d_med+fmt.std_scale*d_std ];
fmt.scale_min = 0.0;
fmt.scale_max = 1.0;
fmt.colormap = fmt.cm_gray;
fmt.fig_pos = fmt.fig_pos_deg_orient;
fnames.save_rel_name = 'degree_of_orientation';
yet_another_plot(fig_no+1,fnames,par,fmt,integ,processing,d,[],ind);
title(fmt.title_end);

subplot(2,2,3);
fmt.title_end = 'orientation and asym. intensity';
fmt.scale = [ d_med-fmt.std_scale*d_std d_med+fmt.std_scale*d_std ];
fmt.scale_min = [];
fmt.scale_max = [];
fmt.colormap = [];
fmt.fig_pos = fmt.fig_pos_comb_asym;
fnames.save_rel_name = 'orientation_asym_int';
yet_another_plot(fig_no+2,fnames,par,fmt,integ,processing,pl_comb_asym,[],ind);
title(fmt.title_end);

subplot(2,2,4);
fmt.title_end = 'orientation, asym. and sym. intensity';
fmt.scale = [ d_med-fmt.std_scale*d_std d_med+fmt.std_scale*d_std ];
fmt.scale_min = [];
fmt.scale_max = [];
fmt.colormap = [];
fmt.fig_pos = fmt.fig_pos_comb_all;
fnames.save_rel_name = 'orientation_asym_sym_int';
yet_another_plot(fig_no+3,fnames,par,fmt,integ,processing,pl_comb_all,[],ind);
title(fmt.title_end);

%save the frames on the video
frame = getframe(movie_figure);
writeVideo(movie_obj1, frame);
fmt.subplot = 0;


%%%%%%%%%%
function [round_value] = round_10(org_value)

pow_of_ten = 10^(floor(log(org_value)/log(10)) -1);
round_value = floor( org_value / pow_of_ten +0.5) * pow_of_ten;

%%%%%%%%%%
function [round_value] = floor_10(org_value)

pow_of_ten = 10^(floor(log(org_value)/log(10)) -1);
round_value = floor( org_value / pow_of_ten ) * pow_of_ten;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [axis_handle] = ...
    yet_another_plot(fig_no,fnames,par,fmt,integ,processing,data_array,cw_rgb,ind)

% to clear the variable
processing.save_figure = [];
% !!! this option seems to be non-functional as fmt.axis_handle is set
% within this function and not returned to the calling function, i.e. the
% field is lost and this function never becomes active. 
% Due to the parfor the plot is not shown anyway. Either one repairs it to
% have the option to use it with for instead of parfor or one removes the
% corresponding non-functional code elements. 
replace_image = ((isfield(fmt,'axis_handle')) && (fmt.axis_handle > 0));
fig_subplot = (isfield(fmt,'subplot')) && (fmt.subplot);

if (~replace_image)
    if (~fig_subplot)
        if isfield(fmt,'figh')
            processing.save_figure = fmt.figh;
        else
            processing.save_figure = figure(fig_no);
            clf;
        end
        hold off;
        % print as layed out on the screen, i.e., preserve aspect ratio
        set(gcf,'PaperPositionMode','auto');
        % paper size
        set(gcf,'PaperType','A4');
        % background color
        set(gcf,'Color','white');
        % resize and position
        set(gcf,'Position',fmt.fig_pos);
    end
    % Font size of axes labels
    set(gca,'FontSize',fmt.axes_font_size);
end

if ~isfield(fmt,'plot_subfigure')
    fmt.plot_subfigure = [];
end


if ((isfield(fmt,'histplot')) && (fmt.histplot))
    % 1D plots like histograms etc.:
    % bins and data are passed in a cell structure as this is parfor
    % compatible
    plot(data_array{1},data_array{2});
    
    % store handle for later use if this is indicated by the existance of
    % the field 
    if (isfield(fmt,'axis_handle_1d'))
        fmt.axis_handle = gca;
    end
else
    % 2D plots:
    if (size(data_array,3) == 3)
        % add scale bar at given point in data points in given length and width
        % in mm
        if (isfield(fmt,'scale_bar_pos_pt'))
            data_array(fmt.scale_bar_pos_pt(2):(fmt.scale_bar_pos_pt(2)+round(fmt.scale_bar_y_mm/par.y_scale)),...
                       fmt.scale_bar_pos_pt(1):(fmt.scale_bar_pos_pt(1)+round(fmt.scale_bar_x_mm/par.x_scale)),:) = 1.0;
        end

        % display 2D array of RGB data, interpolate and orient as specified
        if (fmt.interp_to > 0)
            data_interp_1 = interp2(data_array(:,:,1),fmt.interp_to);
            data_interp = zeros(size(data_interp_1,1),...
                size(data_interp_1,2),3);
            data_interp(:,:,1) = data_interp_1;
            data_interp(:,:,2) = interp2(data_array(:,:,2),fmt.interp_to);
            data_interp(:,:,3) = interp2(data_array(:,:,3),fmt.interp_to);
        else
            data_interp = data_array;
        end
        %if ((isfield(fmt,'fliplr')) && (fmt.fliplr))
        %    data_interp = flipdim(data_interp,2);
        %end
        %if ((isfield(fmt,'flipud')) && (fmt.flipud))
        %    data_interp = flipdim(data_interp,1);
        %end
        if ((isfield(fmt,'permutedim')) && (fmt.permutedim))
            dim_order = [ 2 1 3 ];
            data_interp = permute(data_interp,dim_order);
            x_data = fmt.y;
            y_data = fmt.x;
        else
            x_data = fmt.x;
            y_data = fmt.y;
        end
        if (replace_image)
            % update an existing plot
            set(fmt.image_handle,'CData', data_interp);
        else
            % plot the image
            if ~isfield(fmt,'figh')
                ah = imagesc(data_interp,'XData',x_data,'YData',y_data); %orientation
            else
                fmt.imgh.CData = data_interp;
                fmt.imgh.XData = x_data;
                fmt.imgh.YData = y_data;
                figure(fmt.imgh.Parent.Parent)
                axis auto
            end
            % store handle for later use if this is indicated by the existance 
            % of the field    axis auto
            if (isfield(fmt,'axis_handle'))
                fmt.axis_handle = ah;
            end
        end
    else
        % add scale bar at given point in data points in given length and width
        % in mm
        if (isfield(fmt,'scale_bar_pos_pt'))
            data_array(fmt.scale_bar_pos_pt(2):(fmt.scale_bar_pos_pt(2)+round(fmt.scale_bar_y_mm/par.y_scale)),...
                       fmt.scale_bar_pos_pt(1):(fmt.scale_bar_pos_pt(1)+round(fmt.scale_bar_x_mm/par.x_scale))) = ...
                       2*max(max(data_array));
        end

        % display 2D array of 1D data, orient as specified
        if ((isfield(fmt,'interp_to')) && (fmt.interp_to > 0))
            data_interp = interp2(data_array,fmt.interp_to);
        else
            data_interp = data_array;
        end
        
        if ((isfield(fmt,'permutedim')) && (fmt.permutedim))
            dim_order = [ 2 1 ];
            data_interp = permute(data_interp,dim_order);
            x_data = fmt.y;
            y_data = fmt.x;
        else
            x_data = fmt.x;
            y_data = fmt.y;
        end

        % determine intensity scale
        scale_range = fmt.scale;
        if (~isempty(fmt.scale_min) && (scale_range(1) < fmt.scale_min))
            scale_range(1) = fmt.scale_min;
        end
        if (~isempty(fmt.scale_max) && (scale_range(2) > fmt.scale_max))
            scale_range(2) = fmt.scale_max;
        end
        if (scale_range(1) >= scale_range(2))
            scale_range(2) = scale_range(1) + 0.001;
        end

        if (fig_subplot)
            data_plot = zeros(size(data_interp,1),size(data_interp,2),3);
            data_plot(:,:,1) = (data_interp - scale_range(1)) / (scale_range(2)-scale_range(1));
            data_plot(:,:,2) = data_plot(:,:,1);
            data_plot(:,:,3) = data_plot(:,:,1);
            if (replace_image)
                % update an existing plot
                set(fmt.image_handle,'CData', data_plot);
            else
                ah = imagesc(data_plot,'XData',x_data,'YData',y_data);
                % store handle for later use if this is indicated by the existance 
                % of the field
                if (isfield(fmt,'axis_handle'))
                    fmt.axis_handle = ah;
                end
            end
        else
            if (replace_image)
                % update an existing plot
                set(fmt.image_handle,'CData', data_interp);
            else
                if ~isempty(fmt.plot_subfigure)
                    subplot(fmt.plot_subfigure(1), fmt.plot_subfigure(2), fmt.plot_subfigure(3));
                end
                ah = imagesc(x_data,y_data,data_interp); % symmetric intensity: black white
                % store handle for later use if this is indicated by the existance 
                % of the field
                if (isfield(fmt,'axis_handle'))
                    fmt.axis_handle = ah;
                end
            end

            % set intensity scale
            if (~any(isnan(scale_range)))
                caxis( scale_range );
            end
            
            % display color bar
            colorbar;
        end

    end

    if ((~isfield(fmt,'axis_xy')) || (fmt.axis_xy))
        axis xy;
    else
        axis ij;
    end
    axis equal;
    axis tight;
    if ((~isfield(fmt,'axis_on')) || (fmt.axis_on))
        axis on;
        xlabel('x [ mm ]');
        ylabel('y [ mm ]');
    else
        axis off;
    end
    
    % set return argument: save axes handle of 2D plot for later use
    if (processing.plot_data)
        axis_handle = gca;
    else
        axis_handle = [];
    end
    
end


if ((~isfield(fmt,'title_on')) || (fmt.title_on))
    if ((isfield(fmt,'d_spacing_center_pix')) && (fmt.d_spacing_center_pix))
        d_spacing_pix = mean(par.r_sum{ind});
    else
        d_spacing_pix = fliplr(par.r_sum{ind});
    end
    title_str = [ fmt.title_user_start fmt.title_base ...
        d_spacing_str_par(d_spacing_pix,par,fmt, integ) ];
    if (isfield(fmt,'q_format'))
        title_str = [ title_str ', ' ...
            q_str_par(fliplr(d_spacing_pix),par,fmt) ];
    end
    if (par.r_bgr(ind,1) > 0)
        title_str = [ title_str ' (bgr. sub.)' ];   
    end
    title_str = [ title_str fmt.title_end ];
    title(title_str);
end
if ~isempty(fmt.colormap)
    colormap(fmt.colormap);
end

% resize and position
if ((~isfield(fmt,'subplot')) || (~fmt.subplot))
    set(gcf,'Position',fmt.fig_pos);
end

% display color wheel
if ((~isempty(cw_rgb)) && (~replace_image))
    ap = get(gca,'Position');
    if (ap(3) > 0.04)
        ap(3) = ap(3) - 0.04;
        set(gca,'Position',ap);
    end

    % display a colour wheel outside the plot
    hold on;
    if (isfield(fmt,'cw_size'))
        cw_size = fmt.cw_size;
    else
        cw_size = 0.16;
    end
    if (isfield(fmt,'cw_xoffs'))
        cw_xoffs = fmt.cw_xoffs;
    else
        cw_xoffs = -0.03;
    end
    if (isfield(fmt,'cw_yoffs'))
        cw_yoffs = fmt.cw_yoffs;
    else
        cw_yoffs = -0.03;
    end
    axes('Position', [ap(1) + ap(3) + cw_xoffs - cw_size/2, ap(2) + ap(4) - cw_size + cw_yoffs,...
        cw_size cw_size]);

    imagesc(cw_rgb); %colorwheel plotting
    if ((~isfield(fmt,'axis_xy')) || (fmt.axis_xy))
        axis xy;
    else
        axis ij;
    end
    axis equal;
    axis tight;
    axis off;
end

        
% print and save figure
if (processing.print_fig)
    print('-dpsc',['-P' par.printer_name ]);
end
if (processing.save_fig)
    fmt.file_format = 1;
    dspacing_ftext = ...
        strrep( ['stdscale' num2str(fmt.std_scale,'%4.2f') '_' d_spacing_str_par(fliplr(par.r_sum{ind}),par,fmt,integ) ],' ','');
    fmt.file_format = 0;
    dspacing_ftext = strrep(dspacing_ftext,'.','p');
    dspacing_ftext = [ dspacing_ftext '_' ];
    if (par.r_bgr(ind,1) > 0)
        dspacing_ftext = [ dspacing_ftext 'bgrsub_' ];
    end
    figname = sprintf('%s%s__%s', dspacing_ftext, fnames.save_rel_name, cell2mat(fnames.save_name_suffix(ind)));
    figdir = fullfile(fnames.basedir_figures, '/eps/');
    fig_filename = fullfile( figdir, figname );
    fprintf('saving %s.eps\n',fig_filename);
    print(processing.save_figure, '-depsc','-r1200',[fig_filename '.eps']);
    
    figdir = fullfile( fnames.basedir_figures, '/jpg/');
    fig_filename = fullfile( figdir, figname );
    fprintf('saving %s.jpg\n',fig_filename);
    print(processing.save_figure,'-djpeg','-r300',[fig_filename '.jpg'] );
    
    figdir = fullfile( fnames.basedir_figures, '/fig/');
    fig_filename = fullfile( figdir, figname);
    fprintf('saving %s.fig\n',fig_filename);
    hgsave(processing.save_figure, [fig_filename '.fig'] );
end
if (~processing.plot_data)
    close(gcf);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare the peak fitting or simple peak analysis output data for plotting
function [ fitpeak ] = prepare_peak_fit_plot_data(fit_peak, ind, q, fmt_std_scale)
% peak position:
% initialize
fitpeak.pos = fit_peak.pos(:,:,ind);
% put 2D array in 1D one
fpp_lin = reshape(fitpeak.pos,size(fitpeak.pos,1)*size(fitpeak.pos,2),1);
% convert from pixels to nm, ignore division by zero issues
fpp_lin = 0.1 * 2*pi./interp1(1:length(q),q,fpp_lin);
% set undefined values to zero
fpp_lin(~isfinite(fpp_lin)) = 0;
% put 1D back in 2D array
fitpeak.pos(:) = fpp_lin;
% calculate standard deviation and median for auto-scaling of the
% plot
fitpeak.fpp_std = std(fpp_lin);
fitpeak.fpp_med = median(fpp_lin);

% peak amplitude
% initialize
fitpeak.ampl = fit_peak.ampl(:,:,ind);
% put 2D array in 1D one
fpa_lin = reshape(fitpeak.ampl,size(fitpeak.ampl,1)*size(fitpeak.ampl,2),1);
% set undefined values to zero
fpa_lin(~isfinite(fpa_lin)) = 0;
% put 1D back in 2D array
fitpeak.ampl(:) = fpa_lin;
% calculate standard deviation and median for auto-scaling of the
% plot
fitpeak.fpa_std = std(fpa_lin);
fitpeak.fpa_med = median(fpa_lin);

% peak amplitude:
% initialize
fitpeak.width = fit_peak.width(:,:,ind);
% put 2D array in 1D one
fpw_lin = reshape(fitpeak.width,size(fitpeak.width,1)*size(fitpeak.width,2),1);
% convert from pixels to nm, ignore division by zero issues
fpw_lin = 0.1 * 2*pi./interp1(1:length(q),q,fpw_lin);
% set undefined values to zero
fpw_lin(~isfinite(fpw_lin)) = 0;
% put 1D back in 2D array
fitpeak.width(:) = fpw_lin;
% calculate standard deviation and median for auto-scaling of the
% plot
fitpeak.fpw_std = std(fpw_lin);
fitpeak.fpw_med = median(fpw_lin);


% combined plot all:
fitpeak.pl_comb_all = zeros(size(fit_peak.pos,1),size(fit_peak.pos,2),3);

% scale the hue with the position
fitpeakpos_01 = (fitpeak.pos - fitpeak.fpp_med) / (2 * fmt_std_scale * fitpeak.fpp_std) + 0.5;
fitpeakpos_01( fitpeakpos_01 < 0.0 ) = 0.0;
fitpeakpos_01( fitpeakpos_01 > 1.0 ) = 1.0;
fitpeak.pl_comb_all(:,:,1) = fitpeakpos_01;

% scale the saturation with the peak width
fitpeakwidth_01 = (fitpeak.width - fitpeak.fpw_med) / (2 * fmt_std_scale * fitpeak.fpw_std) + 0.5;
fitpeakwidth_01( fitpeakwidth_01 < 0.0 ) = 0.0;
fitpeakwidth_01( fitpeakwidth_01 > 1.0 ) = 1.0;
fitpeak.pl_comb_all(:,:,2) = fitpeakwidth_01;

% scale the value with the amplitude
% fitpeakampl_01 = (fitpeak.ampl - fitpeak.fpa_med) / (2 * fmt_std_scale * fitpeak.fpa_std) + 0.5;
fitpeakampl_01 = fitpeak.ampl / (fmt_std_scale * fitpeak.fpa_std);
fitpeakampl_01( fitpeakampl_01 < 0.0 ) = 0.0;
fitpeakampl_01( fitpeakampl_01 > 1.0 ) = 1.0;
fitpeak.pl_comb_all(:,:,3) = fitpeakampl_01;

fitpeak.pl_comb_all = hsv2rgb(fitpeak.pl_comb_all);


    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the peak fitting or simple peak analysis outputdata
function [ ah ] = plot_peak_fit_data(fig_peak_pos,fig_peak_ampl,fig_peak_width, fig_peak_all, ...
    fmt_fig_pos_fit_peak, ...
    fnames,par,fmt,integ,processing, ...
    basename, ...
    fitpeak, ...
    ind)

% axes handles
ah = zeros(4,1);

if fig_peak_pos
    fmt.title_end = sprintf(', %s, 2\\pi/(peak pos.) [nm]',basename);
    fmt.scale = [ fitpeak.fpp_med-fmt.std_scale*fitpeak.fpp_std fitpeak.fpp_med+fmt.std_scale*fitpeak.fpp_std ];
    fmt.scale_min = 0;
    fmt.scale_max = [];
    fmt.colormap = fmt.cm_gray;
    fmt.fig_pos = fmt_fig_pos_fit_peak.pos;
    fnames.save_rel_name = [ basename '_pos' ];
    ah(1) = yet_another_plot(ind*10+110,fnames,par,fmt,integ,processing,fitpeak.pos,[],ind);
end
%%%%% peak amplitude
if fig_peak_ampl
    fmt.title_end = sprintf(', %s, peak ampl. [a.u.]', basename);
    fmt.scale = [ fitpeak.fpa_med-fmt.std_scale*fitpeak.fpa_std fitpeak.fpa_med+fmt.std_scale*fitpeak.fpa_std ];
    fmt.scale_min = 0;
    fmt.scale_max = [];
    fmt.colormap = fmt.cm_gray;
    fmt.fig_pos = fmt_fig_pos_fit_peak.ampl;
    fnames.save_rel_name = [ basename '_ampl' ];
    ah(2) = yet_another_plot(ind*10+111,fnames,par,fmt,integ,processing,fitpeak.ampl,[],ind);
end
%%%%% peak width
if fig_peak_width
    fmt.title_end = sprintf(', %s, 2\\pi/(peak width) [nm]', basename);
    fmt.scale = [ fitpeak.fpw_med-fmt.std_scale*fitpeak.fpw_std fitpeak.fpw_med+fmt.std_scale*fitpeak.fpw_std ];
    fmt.scale_min = 0;
    fmt.scale_max = [];
    fmt.colormap = fmt.cm_gray;
    fmt.fig_pos = fmt_fig_pos_fit_peak.width;
    fnames.save_rel_name = [ basename '_width'];
    ah(3) = yet_another_plot(ind*10+112,fnames,par,fmt,integ,processing,fitpeak.width,[],ind);
end
%%%%% peak amplitude, position, width combined
if fig_peak_all
    fmt.title_end = sprintf(', %s, peak ampl., 2\\pi/pos., 2\\pi/width', basename);
    fmt.scale = [];
    fmt.scale_min = 0;
    fmt.scale_max = [];
    fmt.colormap = [];
    fmt.fig_pos = fmt_fig_pos_fit_peak.all;
    fnames.save_rel_name = [ basename '_ampl_pos_width' ];
    ah(4) = yet_another_plot(ind*10+113,fnames,par,fmt,integ,processing,fitpeak.pl_comb_all,[],ind);
end
