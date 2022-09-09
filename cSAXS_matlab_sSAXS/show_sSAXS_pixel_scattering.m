% function show_sSAXS_pixel_scattering((filename,[[,<name>,<value>] ...])
% Shows the scattering pattern and/or azimuthal integration plots for a
% selected point in the sSAXS colorful pictures. Will only work for fourier
% components computed with code newer than Dec 13, 2017.
%
% Inputs
%    filename       Fourier components filename, set = [] for using GUI open. E.g. '~/Data10/analysis/fourier_components/my_cute_sample_here/sample_id/20026_fourier.mat';
%
% The optional <name>,<value> pairs are:
% 
% 'which_qrange_to_plot'    Which of the possibly multiple processed q-ranges to be shown (default = 1).
% 'fig_type'                Which figure to show ( default = 'fig_asym_sym_int'). 
%                           Options: fig_asym_sym_int, 'fig_asym_int',
%                           'fig_orientation', fig_orient_degree,
%                           fig_orient_degree_histogram, fig_scatt_asym,
%                           fig_cos_dev, fig_I_pow, 
%                           fig_simple_peak.pos, fig_simple_peak.ampl, fig_simple_peak.width, fig_simple_peak.all, 
%                           fig_fit_peak.pos, fig_fit_peak.ampl, fig_fit_peak.width, fig_fit_peak.all
% 'coordinates' = [];       [x y] spatial coordinates on the image, or empty for GUI input
% 'show_diffraction_pattern'    Show the diffraction pattern (default = true)
% 'show_integrated_plot'        Show the plot of the azimuthally integrated data (default = true)
% 'compile_x12sa_filename_args' Extra <name, value> arguments to be passed to utils.compile_x12sa_filename 
%                               to compile the diffraction pattern filename (default = {}). E.g. {'BaseName','e16656_'}
% 'image_show_args'             Extra <name, value> arguments to the function plotting.image_show (default = {}). E.g. {'RowFrom',<0-max>}
% 'plot_radial_integ_args'      Extra <name, value> arguments to the function plotting.plot_radial_integ (default = {}). E.g. {'SegAvg',<0-no, 1-yes>}
% 'basedir_integrated_data'     Directory of the azimuthally integrated data (default = [ ]). If empty this path is taken from fnames.basedir_integrated_data, e.g. '~/Data10/analysis/radial_integration/'
%
% Example: show_sSAXS_pixel_scattering([])

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function show_sSAXS_pixel_scattering(filename,varargin)

% Defaults
%filename = '~/Data10/analysis/fourier_components/test_section35um/mouse2_section50um_onecorner_1/00277_fourier.mat';
% which_qrange_to_plot = 1;   % Which of the possibly multiple processed q-ranges to be shown.
% fig_type = 'fig_asym_sym_int';  % Options: fig_asym_sym_int, 'fig_asym_int', 'fig_orientation', fig_orient_degree, fig_orient_degree_histogram, fig_scatt_asym, fig_cos_dev, fig_I_pow
% coordinates = [];   % [X Y] spatial coordinates on the image, or empty for GUI input
% show_diffraction_pattern = true;
% show_integrated_plot = true;
% compile_x12sa_filename_args = {}; % Arguments to the function utils.compile_x12sa_filename to find the name of the scatterin pattern frame
% image_show_args = {}; % Arguments to the function plotting.image_show 
% plot_radial_integ_args = {}; % Arguments to the function plotting.plot_radial_integ 
% basedir_integrated_data = []; % If empty it is taken from fnames.basedir_integrated_data, e.g. '~/Data10/analysis/radial_integration/'
%%%%%%%%%%%%%

%%% Parse inputs %%%
par = inputParser;
par.addParameter('which_qrange_to_plot',1)
par.addParameter('fig_type', 'fig_asym_sym_int')
par.addParameter('coordinates',[])
par.addParameter('show_diffraction_pattern',true)
par.addParameter('show_integrated_plot',true)
par.addParameter('compile_x12sa_filename_args',{})
par.addParameter('image_show_args',{})
par.addParameter('plot_radial_integ_args',{})
par.addParameter('basedir_integrated_data',[])

par.parse(varargin{:})
r = par.Results;
which_qrange_to_plot = r.which_qrange_to_plot;
fig_type = r.fig_type;
initial_coordinates = r.coordinates;
show_diffraction_pattern = r.show_diffraction_pattern;
show_integrated_plot = r.show_integrated_plot;
compile_x12sa_filename_args = r.compile_x12sa_filename_args; 
image_show_args = r.image_show_args; 
initial_plot_radial_integ_args = r.plot_radial_integ_args;
basedir_integrated_data = r.basedir_integrated_data; 

% addpath ../   % Relies on matlab base package
if isempty(filename)
    [FileName,PathName,FilterIndex] = uigetfile('','Select Fourier components file');
    filename = fullfile(PathName,FileName);
end
fprintf('Loading %s\n',filename);
fprintf('Click on the sample image on the pixel you want to examine.\n');
load(filename)

% Reset figure types to zero, avoid plotting several figures
four_dat.fmt.fig_sym_int                     = 0;
four_dat.fmt.fig_asym_int                    = 0;
four_dat.fmt.fig_orientation                 = 0;
four_dat.fmt.fig_orient_degree               = 0;
four_dat.fmt.fig_orient_degree_histogram     = 0;
four_dat.fmt.fig_asym_sym_int                = 0;
four_dat.fmt.fig_scatt_asym                  = 0;
four_dat.fmt.fig_cos_dev                     = 0;
four_dat.fmt.fig_I_pow                       = 0;
% simple peak analysis
four_dat.fmt.fig_simple_peak.pos = 0;       % peak position
four_dat.fmt.fig_simple_peak.ampl = 0;      % peak amplitude
four_dat.fmt.fig_simple_peak.width = 0;     % peak width
four_dat.fmt.fig_simple_peak.all = 0;       % peak amplitude, position, width
% peak fitting
four_dat.fmt.fig_fit_peak.pos = 0;       % peak position
four_dat.fmt.fig_fit_peak.ampl = 0;      % peak amplitude
four_dat.fmt.fig_fit_peak.width = 0;     % peak width
four_dat.fmt.fig_fit_peak.all = 0;       % peak amplitude, position, width

four_dat.fmt.parallel_plotting = false; % Deactivate plotting with parfor so that the figure will actually show

% activate the user preferred plot (trusting that only valid field-names
% are specified as parameter)
field_names = split(fig_type,'.');
if (length(field_names) == 1)
    four_dat.fmt.(fig_type) = 1;
else
    four_dat.fmt.(field_names{1}).(field_names{2}) = 1;
end
% No interpolation for this figure
four_dat.fmt.interpol = 0;

% Reset processing properties, avoid resaving figures and jpegs
processing.load_data    = 0;
processing.save_fig     = 0;
processing.plot_data    = 1;
processing.movie_fig    = 0;
processing.print_fig    = 0;
processing.which_qrange_to_plot = which_qrange_to_plot; % Set to empty or if it does not exist then it plots all


% If coordinates have not been specified as command line argument then
% repeatedly wait for a mouse-click specifying the location to display
% information on. 
clear_fig = 1;
coordinates = [];
while (1)
    % initially and upon clear-figure request redraw the Fourier result
    % figure
    if (clear_fig)
        plot_fourier_result(four_dat, processing);
        fig_handle = gcf;
        % in case of a clear-figure request restore the lines marking the
        % latest coordinates
        if (~isempty(coordinates))
            plot_lines(fig_handle,coordinates);
        end
    end
    
    % coordinates are either mouse-click or command-line parameter
    % specified
    if isempty(initial_coordinates)
        fprintf('Left button -- add I(q) plot, center button -- clear figure for I(q) plot, right button -- exit\n');
        figure(fig_handle);
        [coordinates(1), coordinates(2), mouse_button] = ginput(1);
        % right button: exit
        if (mouse_button > 2)
            break;
        end
        % left-click: keep I(q) plot, middle-click: clear plot 
        if (mouse_button == 1)
            clear_fig = 0;
        else
            clear_fig = 1;
        end
    else
        coordinates = initial_coordinates;
    end
    % convert (x,y) coordinate back to integer indices to the data arrays
    indices(1) = round(coordinates(1)/four_dat.par.x_scale);
    indices(2) = round(coordinates(2)/four_dat.par.y_scale);
    indices = indices+1;
    % check for clicks outside the plot (axes area)
    clip_indices = max(indices,[1 1]);
    clip_indices = min(clip_indices, fliplr(size(four_dat.scan_num)));
    if (any(clip_indices ~= indices))
        fprintf('Coordinates outside figure -- relocating into the region confined by the axes.\n');
        indices = clip_indices;
        coordinates(1) = (indices(1)-1) * four_dat.par.x_scale;
        coordinates(2) = (indices(2)-1) * four_dat.par.y_scale;
    end
    % print coordinates and mark position with two lines, i.e. a cross
    fprintf('(x,y)       = (%.3f, %.3f) mm     ',coordinates)
    fprintf('(indx,indy) = (%d, %d)\n',indices)
    plot_lines(fig_handle,coordinates);
    
    % get scan-line number and point within this line number
    scan_num    = four_dat.scan_num(indices(2),indices(1));
    scan_point  = four_dat.scan_point(indices(2),indices(1));

    if show_diffraction_pattern
        % Try first the most common, a subexposure in a continuous line
        num_arg_comp = numel(compile_x12sa_filename_args);
        compile_x12sa_filename_args{num_arg_comp+1} = 'SubExpNo';
        compile_x12sa_filename_args{num_arg_comp+2} = scan_point;
        image_filename = utils.compile_x12sa_filename(scan_num,0,compile_x12sa_filename_args);
        if ~exist(image_filename,'file')
            % Try to find the index as scan point
            warning('Did not find %s, trying as a point index',image_filename);
            image_filename = utils.compile_x12sa_filename(scan_num,scan_point,compile_x12sa_filename_args);
        end
        plotting.image_show(image_filename,image_show_args);
    end

    if show_integrated_plot
        if isempty(basedir_integrated_data)
            basedir_integrated_data = four_dat.fnames.basedir_integrated_data;
        end
        plot_radial_integ_args = initial_plot_radial_integ_args;
        num_radplot_args = numel(plot_radial_integ_args);
        plot_radial_integ_args{num_radplot_args+1} = 'PointRange';
        plot_radial_integ_args{num_radplot_args+2} = scan_point+1;
        % pass on clear-figure request 
        plot_radial_integ_args{num_radplot_args+3} = 'ClearFig';
        if (clear_fig)
            plot_radial_integ_args{num_radplot_args+4} = 1;
        else
            plot_radial_integ_args{num_radplot_args+4} = 0;
        end
        radial_integ_filename = sprintf(fullfile(basedir_integrated_data,four_dat.fnames.basename_integrated_data),scan_num);
        % plot I(q) or pixel
        plotting.plot_radial_integ(radial_integ_filename,plot_radial_integ_args);
    end
    
    if (~isempty(initial_coordinates))
        break;
    end
    fprintf('\n')
end    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_lines(fig_handle,coordinates)
% find axes handle that points to an object with non-empty title string,
% i.e. find the axes that belong to the plot(rather than for example the
% color-disc legend)
all_axes = findobj(fig_handle);
if (isempty(all_axes))
    fprintf('No axes handles found.\n');
    return;
end
axes_handle = [];
for ind = 1:length(all_axes)
    if ((isprop(all_axes(ind),'Title')) && (isprop(all_axes(ind).Title,'String')) && (~isempty(all_axes(ind).Title.String)))
        axes_handle = all_axes(ind);
        break;
    end
end
if (isempty(axes_handle))
    fprintf('No axis handle with non-empty title string found.\n');
end
% finally make the plot the current axes
axes(axes_handle);
% and ensure that the lines do not replace the current plot
hold(axes_handle,'on');

% plot vertical and horizontal line
% plotting.vline(coordinates(1),'g');
% plotting.hline(coordinates(2),'g');

% vline and hline work fine with a limited set of parameters
% to access other parameters like arbitrary line-colors it is easiest to
% plot the line here directly
ylim=get(gca,'ylim');
plot([coordinates(1) coordinates(1)],ylim,'Color',[0.5 0.5 0.5]);
xlim=get(gca,'xlim');
plot(xlim, [coordinates(2) coordinates(2)],'Color',[0.5 0.5 0.5]);
