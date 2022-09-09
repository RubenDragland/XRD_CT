% comparison of the result from optimization and measured data
% compare the 2D projections

close all
clear all

addpath ../scanning_saxs/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_path = '/das/work/units/csaxs/p17283/Bone_sample/';%'~/Data10';         % = ~/Data10 for online analysis, provide the path for offline analysis Ex: '/das/work/p16/p16649/'
sample_name = 'bone';%'sample_name';    % name given in the saxs_caller_template
add_name = 'test_debug_20190627';%'ID';         % additional name the optimizations: = [ ] if not needed
include_IRTT = 0;               % = 1 to include the results from IRTT
                                % = 0 to not include IRTT    
theta_det = pi/2;               % Ewald sphere correction, for SAXS pi/2
which_projections = [1:10:200]; % which projections to plot, example [1:10:300]
plot_corrected_data = 1;        % =0, if 2D projections don't have corrections
                                % =1, if corrections have been applied during alignment
apply_mask = 1;                 % = 0 does not apply mask
                                % = 1 apply mask
std_scale = [];                 % = [] does not change  the std scale
                                % = 3, changes to scale to 3*std of the data (as in SAXS_caller_template)
cw_size = [];                   % = [], color wheel will be same size (about 0.16)
                                % = 0.09, for decreasing  the size of  the
                                % color wheel
save_video = 0;                 % = 1, to save video, = 0 does not save
save_figure = 0;                % = 1, to save figure, = 0 does not save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load aligned projections
filename = fullfile(base_path,sprintf('analysis/SASTT/%s/projection_data/SASTT_%s_aligned_ASTRA.mat', ...
    sample_name, sample_name));
load(filename);


% prepare video
if save_video
    videoname = sprintf('video_%s_%s.avi', sample_name, add_name);
    videoname = fullfile(base_path,  sprintf('analysis/SASTT/%s/SH/%s/figures/%s', sample_name, add_name, videoname));
    if ~(exist('videoname', 'file'))
        v = VideoWriter(videoname);
        v.FrameRate = 2;
        open(v)
    end
end
screen_size = get(0,'screensize');
four_dat.fmt.fig_pos_sym_int =  [150 150 screen_size(3)-300 screen_size(4)-300 ];
four_dat.fmt.fig_pos_comb_all =  [150 150 screen_size(3)-300 screen_size(4)-300 ];

%% loop over THE projections
for kk = 1:numel(which_projections)
    % select one projection
        one_projection = projection(which_projections(kk));
        one_projection.integ.theta_det = theta_det;
        %comparing
        fprintf('Comparing projection %05d at rotX %.1f and rotY %.1f \n', which_projections(kk), one_projection.rot_x, one_projection.rot_y);
        fprintf('** \n');
        
    if ~(include_IRTT)
        % plot the original data
        % load the measured 2D projections
        filename = sprintf('%sanalysis/fourier_components/%s/%s%05d_fourier.mat', base_path, ...
            sample_name, one_projection.fnames.rel_dir, one_projection.fnames.first_scan_no);
        load(filename)
        
        % if needed, change  the saturation of the images
        if ~(isempty(std_scale))
            four_dat.fmt.std_scale = std_scale;
        end
        % if needed, change  the size of the color wheel
        if ~(isempty(cw_size))
            four_dat.fmt.cw_size = cw_size;
        end
        
        % needed to only plot the data
        processing.load_data = 0; % 1 to force re-load and re-process
        processing.save_fig  = 0;
        processing.plot_data = 1;
        processing.print_fig = 0;
        processing.movie_fig = 0;
        % define which q-range to plot
        processing.which_qrange_to_plot = four_dat.par.which_qrange;
        
        % only plot the important plots
        four_dat.fmt.fig_sym_int = 1;        % symmetric intensity over all segments: f1_amp
        four_dat.fmt.fig_asym_sym_int = 0;   % combination of orientation (colorwheel), degree of orientation (hue)
        four_dat.fmt.fig_asym_int = 0;       % amplitude of the anisotropy: f2_amp
        four_dat.fmt.fig_orientation = 0;    % angle of orientation (in degrees): f2_phase
        four_dat.fmt.fig_orient_degree = 0;  % degree of orientation from f2_amp/f1_amp
        four_dat.fmt.fig_orient_degree_histogram = 0; % histogram of the degree of orientation vs. orientation
        four_dat.fmt.fig_scatt_asym = 0;     % scattering asymmetry
        four_dat.fmt.fig_cos_dev = 0;        % how the fit of azimuthal integration deviates from a cosine
        four_dat.fmt.fig_I_pow = 0;          % slopes of I(q) in the selected q-region
        four_dat.fmt.fig_peak_pos = 0;       % peak position and amplitude
        % define positions

        
        % turn title off
        four_dat.fmt.title_on = 0;
        
        
        if plot_corrected_data
            processing.which_qrange_to_plot = 1;
            % calculate the fourier coefficients
            [four_dat] = calculate_four_dat(four_dat, one_projection.data);
        end
        
        % if we need to apply the mask
        if (apply_mask)
            if ~(plot_corrected_data)
                % adjust sizes, if mask and original data have different
                % sizes
                temp = [];
                names = fieldnames(four_dat);
                for pp = 2:9 % only the fields needed here
                    % y dimension
                    while size(four_dat.(names{pp}), 1) < size(one_projection.window_mask, 1)
                        temp = [four_dat.(names{pp})(1, :, :); four_dat.(names{pp})];
                        four_dat.(names{pp}) = temp;
                    end
                    % x dimension
                    while size(four_dat.(names{pp}), 2) < size(one_projection.window_mask, 2)
                        temp = [four_dat.(names{pp})(:, 1, :), four_dat.(names{pp})];
                        four_dat.(names{pp}) = temp;
                    end
                end
            end
            % compile return structure for plotting and mask
            four_dat.f1_amp = four_dat.f1_amp.*one_projection.window_mask;
            four_dat.f2_amp = four_dat.f2_amp.*one_projection.window_mask;
            four_dat.f2_phase = four_dat.f2_phase.*one_projection.window_mask;
            four_dat.orientation = four_dat.orientation.*one_projection.window_mask;
            four_dat.degree_orient = four_dat.degree_orient.*one_projection.window_mask;
        end
        
        four_dat.fmt.parallel_plotting = false;
        
        % plot the mesured symmetric intensity
        figh = figure(100);
        clf(figh)
        axesh = gca;
        subplot(2,2,1)
        imgh = imagesc(1);
        four_dat.fmt.figh = figh;
        four_dat.fmt.imgh = imgh;
        % plot
        plot_fourier_result(four_dat, processing);
        title(sprintf('Measured Sym. Int., rotX %.1f, rotY %.1f', one_projection.rot_x, one_projection.rot_y), 'FontSize', 15);
        
        
        % plot the measured orientation
        four_dat.fmt.fig_sym_int = 0;        % symmetric intensity over all segments: f1_amp
        four_dat.fmt.fig_asym_sym_int = 1;   % combination of orientation (colorwheel), degree of orientation (hue)
        figh = figure(100);
        subplot(2,2,3)
        imgh = imagesc(1);
        axesh = gca;
        four_dat.fmt.figh = figh;
        four_dat.fmt.imgh = imgh;
        % plot
        plot_fourier_result(four_dat, processing);
        title(axesh, sprintf('Measured Orient., rotX %.1f, rotY %.1f', one_projection.rot_x, one_projection.rot_y), 'FontSize', 15);
        
        % Calculate 2D projections from the optimized data with SH
        filename = fullfile(base_path, sprintf('analysis/SASTT/%s/SH/%s/results/', sample_name, add_name));
        filename = fullfile(filename, sprintf('result_%s_q%d-%d_all_again_%s.mat', sample_name, projection(1).par.r_sum{1}(1), ...
            projection(1).par.r_sum{1}(end), add_name));
        load(filename)
        
        % For more efficiency the loop over projections should be done
        % inside SAXS_tomo_3D_err_metric
        [~,~,proj_SH] = optimization.SAXS_tomo_3D_err_metric([],p, s, one_projection);
        % load data as the variables must be reprocessed
        processing.load_data = 1; % 1 to force re-load and re-process
        processing.plot_data = 1;
        % define which q-range to plot
        processing.which_qrange_to_plot = 1;
        
        % calculate the fourier coefficients
        [four_dat] = calculate_four_dat(four_dat, proj_SH.projection);
        
        % if needed to apply the mask
        if (apply_mask)
            % compile return structure for plotting and mask
            four_dat.f1_amp = four_dat.f1_amp.*one_projection.window_mask;
            four_dat.f2_amp = four_dat.f2_amp.*one_projection.window_mask;
            four_dat.f2_phase = four_dat.f2_phase.*one_projection.window_mask;
            four_dat.orientation = four_dat.orientation.*one_projection.window_mask;
            four_dat.degree_orient = four_dat.degree_orient.*one_projection.window_mask;
        end
        
        % plot the simulated symmetric intensity
        four_dat.fmt.fig_sym_int = 1;        % symmetric intensity over all segments: f1_amp
        four_dat.fmt.fig_asym_sym_int = 0;   % combination of orientation (colorwheel), degree of orientation (hue)
        figh = figure(100);
        subplot(2,2,2)
        imgh = imagesc(1);
        four_dat.fmt.figh = figh;
        four_dat.fmt.imgh = imgh;
        % plot
        plot_fourier_result(four_dat, processing);
        title(sprintf('Simulated Sym. Int., rotX %.1f, rotY %.1f', one_projection.rot_x, one_projection.rot_y) , 'FontSize', 15);
        
        % plot the simulated orientation
        four_dat.fmt.fig_sym_int = 0;        % symmetric intensity over all segments: f1_amp
        four_dat.fmt.fig_asym_sym_int = 1;   % combination of orientation (colorwheel), degree of orientation (hue)
        figh = figure(100);
        subplot(2,2,4)
        imgh = imagesc(1);
        axesh = gca;
        
        four_dat.fmt.figh = figh;
        four_dat.fmt.imgh = imgh;
        % plot
        plot_fourier_result(four_dat, processing);
        title(axesh, sprintf('Simulated Orient., rotX %.1f, rotY %.1f', one_projection.rot_x, one_projection.rot_y), 'FontSize', 15);
        
    else
        %% plot the original data
        % load the measured 2D projections
        filename = sprintf('%sanalysis/fourier_components/%s/%s%05d_fourier.mat', base_path, ...
            sample_name, one_projection.fnames.rel_dir, one_projection.fnames.first_scan_no);
        load(filename)
        
        % if needed, change  the saturation of the images
        if ~(isempty(std_scale))
            four_dat.fmt.std_scale = std_scale;
        end
        % if needed, change  the size of the color wheel
        if ~(isempty(cw_size))
            four_dat.fmt.cw_size = cw_size;
        end
        
        % needed to only plot the data
        processing.load_data = 0; % 1 to force re-load and re-process
        processing.save_fig  = 0;
        processing.plot_data = 1;
        processing.print_fig = 0;
        processing.movie_fig = 0;
        % define which q-range to plot
        processing.which_qrange_to_plot = four_dat.par.which_qrange;
        
        % only plot the important plots
        four_dat.fmt.fig_sym_int = 1;        % symmetric intensity over all segments: f1_amp
        four_dat.fmt.fig_asym_sym_int = 0;   % combination of orientation (colorwheel), degree of orientation (hue)
        four_dat.fmt.fig_asym_int = 0;       % amplitude of the anisotropy: f2_amp
        four_dat.fmt.fig_orientation = 0;    % angle of orientation (in degrees): f2_phase
        four_dat.fmt.fig_orient_degree = 0;  % degree of orientation from f2_amp/f1_amp
        four_dat.fmt.fig_orient_degree_histogram = 0; % histogram of the degree of orientation vs. orientation
        four_dat.fmt.fig_scatt_asym = 0;     % scattering asymmetry
        four_dat.fmt.fig_cos_dev = 0;        % how the fit of azimuthal integration deviates from a cosine
        four_dat.fmt.fig_I_pow = 0;          % slopes of I(q) in the selected q-region
        four_dat.fmt.fig_peak_pos = 0;       % peak position and amplitude
        % define positions
        screen_size = get(0,'screensize');
%         four_dat.fmt.fig_pos_sym_int =  [150 150 screen_size(3)-300 screen_size(4)-300 ];
%         four_dat.fmt.fig_pos_comb_all =  [150 150 screen_size(3)-300 screen_size(4)-300 ];
        % turn title off
        four_dat.fmt.title_on = 0;
        
        
        % if SAXS data has been corrected
        if plot_corrected_data
            processing.which_qrange_to_plot = 1;
            % calculate the fourier coefficients
            [four_dat] = calculate_four_dat(four_dat, one_projection.data);
        end
        
        % if needed to apply a mask
        if (apply_mask)
            if ~(plot_corrected_data)
                % adjust sizes, if mask and original data have different
                % sizes
                temp = [];
                names = fieldnames(four_dat);
                for pp = 2:9 % only the fields needed here
                    % y dimension
                    while size(four_dat.(names{pp}), 1) < size(one_projection.window_mask, 1)
                        temp = [four_dat.(names{pp})(1, :, :); four_dat.(names{pp})];
                        four_dat.(names{pp}) = temp;
                    end
                    % x dimension
                    while size(four_dat.(names{pp}), 2) < size(one_projection.window_mask, 2)
                        temp = [four_dat.(names{pp})(:, 1, :), four_dat.(names{pp})];
                        four_dat.(names{pp}) = temp;
                    end
                end
            end
            % compile return structure for plotting and mask
            four_dat.f1_amp = four_dat.f1_amp.*one_projection.window_mask;
            four_dat.f2_amp = four_dat.f2_amp.*one_projection.window_mask;
            four_dat.f2_phase = four_dat.f2_phase.*one_projection.window_mask;
            four_dat.orientation = four_dat.orientation.*one_projection.window_mask;
            four_dat.degree_orient = four_dat.degree_orient.*one_projection.window_mask;
        end
        
        
        % plot the mesured symmetric intensity
        figh = figure(100);
        clf(figh)
        axesh = gca;
        subplot(2,3,1)
        imgh = imagesc(1);
        four_dat.fmt.figh = figh;
        four_dat.fmt.imgh = imgh;
        % plot
        plot_fourier_result(four_dat, processing);
        title(sprintf('Measured, rotX %.1f, rotY %.1f', one_projection.rot_x, one_projection.rot_y), 'FontSize', 10);
        
        
        % plot the measured orientation
        four_dat.fmt.fig_sym_int = 0;        % symmetric intensity over all segments: f1_amp
        four_dat.fmt.fig_asym_sym_int = 1;   % combination of orientation (colorwheel), degree of orientation (hue)
        figh = figure(100);
        subplot(2,3,4)
        imgh = imagesc(1);
        axesh = gca;
        four_dat.fmt.figh = figh;
        four_dat.fmt.imgh = imgh;
        % plot
        plot_fourier_result(four_dat, processing);
        title(axesh, sprintf('Measured, rotX %.1f, rotY %.1f', one_projection.rot_x, one_projection.rot_y), 'FontSize', 10);
  
        
        % Calculate 2D projections from the optimized data with SH
        filename = fullfile(base_path, sprintf('analysis/SASTT/%s/SH/%s/results/', sample_name, add_name));
        filename = fullfile(filename, sprintf('result_%s_q%d-%d_all_again_%s.mat', sample_name, projection(1).par.r_sum{1}(1), ...
            projection(1).par.r_sum{1}(end), add_name));
        load(filename)
        
        %%%%%%%%%%%%%%%%%%
        % For more efficiency the loop over projections should be done
        % inside SAXS_tomo_3D_err_metric
        [~,~,proj_SH] = SAXS_tomo_3D_err_metric([],p, s, one_projection);
        %%%%%%%%%%%%%%%%%%
        
        % load data as the variables must be reprocessed
        processing.load_data = 1; % 1 to force re-load and re-process
        processing.plot_data = 1;
        % define which q-range to plot
        processing.which_qrange_to_plot = 1;
        
        % calculate the fourier coefficients
        [four_dat] = calculate_four_dat(four_dat, proj_SH.projection);
        
        if (apply_mask)
            % compile return structure for plotting and mask
            four_dat.f1_amp = four_dat.f1_amp.*one_projection.window_mask;
            four_dat.f2_amp = four_dat.f2_amp.*one_projection.window_mask;
            four_dat.f2_phase = four_dat.f2_phase.*one_projection.window_mask;
            four_dat.orientation = four_dat.orientation.*one_projection.window_mask;
            four_dat.degree_orient = four_dat.degree_orient.*one_projection.window_mask;
        end
        
        
        % plot the simulated symmetric intensity
        four_dat.fmt.fig_sym_int = 1;        % symmetric intensity over all segments: f1_amp
        four_dat.fmt.fig_asym_sym_int = 0;   % combination of orientation (colorwheel), degree of orientation (hue)
        figh = figure(100);
        subplot(2,3,2)
        imgh = imagesc(1);
        four_dat.fmt.figh = figh;
        four_dat.fmt.imgh = imgh;
        % plot
        plot_fourier_result(four_dat, processing);
        title(sprintf('SASTT Sym. Int., rotX %.1f, rotY %.1f', one_projection.rot_x, one_projection.rot_y) , 'FontSize', 10);
        
        % plot the simulated orientation
        four_dat.fmt.fig_sym_int = 0;        % symmetric intensity over all segments: f1_amp
        four_dat.fmt.fig_asym_sym_int = 1;   % combination of orientation (colorwheel), degree of orientation (hue)
        figh = figure(100);
        subplot(2,3,5)
        imgh = imagesc(1);
        axesh = gca;
        
        four_dat.fmt.figh = figh;
        four_dat.fmt.imgh = imgh;
        % plot
        plot_fourier_result(four_dat, processing);
        title(axesh, sprintf('SASTT Orient., rotX %.1f, rotY %.1f', one_projection.rot_x, one_projection.rot_y), 'FontSize', 10);
        
        % include IRTT data
        filename = fullfile(base_path, sprintf('analysis/SASTT/%s/IRTT/%s/results/', sample_name, add_name));
        filename = fullfile(filename, sprintf('IRTT_%s_q_%d-%d_%s.mat', sample_name, projection(1).par.r_sum{1}(1), ...
            projection(1).par.r_sum{1}(end), add_name));
        load(filename)
        
        B = squeeze(B_segs(which_projections(kk),:,:));
        % projsize is [segments Y X]
        proj_IR = proj_simul_tensor(tomotensor, one_projection.Rot_exp, B);
        % change dimensions of proj_IR to [Y X segments]
        proj_IR = permute(proj_IR, [2 3 1]);
        % calculate  the fourier data
        [four_dat] = calculate_four_dat(four_dat, proj_IR);
        
        if (apply_mask)
            % compile return structure for plotting and mask
            four_dat.f1_amp = four_dat.f1_amp.*one_projection.window_mask;
            four_dat.f2_amp = four_dat.f2_amp.*one_projection.window_mask;
            four_dat.f2_phase = four_dat.f2_phase.*one_projection.window_mask;
            four_dat.orientation = four_dat.orientation.*one_projection.window_mask;
            four_dat.degree_orient = four_dat.degree_orient.*one_projection.window_mask;
        end
        
        % load data as the variables must be reprocessed
        processing.load_data = 1; % 1 to force re-load and re-process
        processing.plot_data = 1;
        % define which q-range to plot
        processing.which_qrange_to_plot = 1;

        % plot the simulated symmetric intensity
        four_dat.fmt.fig_sym_int = 1;        % symmetric intensity over all segments: f1_amp
        four_dat.fmt.fig_asym_sym_int = 0;   % combination of orientation (colorwheel), degree of orientation (hue)
        figh = figure(100);
        subplot(2,3,3)
        imgh = imagesc(1);
        four_dat.fmt.figh = figh;
        four_dat.fmt.imgh = imgh;
        % plot
        plot_fourier_result(four_dat, processing);
        title(sprintf('IRTT Sym. Int., rotX %.1f, rotY %.1f', one_projection.rot_x, one_projection.rot_y) , 'FontSize', 10);
        
        % plot the simulated orientation
        four_dat.fmt.fig_sym_int = 0;        % symmetric intensity over all segments: f1_amp
        four_dat.fmt.fig_asym_sym_int = 1;   % combination of orientation (colorwheel), degree of orientation (hue)
        figh = figure(100);
        subplot(2,3,6)
        imgh = imagesc(1);
        axesh = gca;
        
        four_dat.fmt.figh = figh;
        four_dat.fmt.imgh = imgh;
        % plot
        plot_fourier_result(four_dat, processing);
        title(axesh, sprintf('IRTT Orient., rotX %.1f, rotY %.1f', one_projection.rot_x, one_projection.rot_y), 'FontSize', 10);
    end
    
    % for saving the results
    if save_figure
        % save the figure
        figname = sprintf('scan_%05d_rotx%.0f_roty%.0f', four_dat.fnames.first_scan_no, one_projection.rot_x, one_projection.rot_y);
        %make a new folder to save  the plotting comparison
        basedir_figures = fullfile(base_path,  sprintf('analysis/SASTT/%s/SH/%s/figures/%s', sample_name, add_name, figname));
        print(gcf, basedir_figures,'-dpng','-r300');
    end
    if save_video
        frame_plot = getframe(gcf);
        writeVideo(v, frame_plot);
    end
     pause(1)
end




%*-------------------------------------------------------------------------------------*
%|                                                                                     |
%|  Except where otherwise noted, this work is licensed under a                        |     
%|                                                                                     |   
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0                          |
%|  International (CC BY-NC-SA 4.0) license.                                           |
%|                                                                                     |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)                  |
%|                                                                                     |
%|      Author: CXS group, PSI                                                         |
%*------------------------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If this code, or subfunctions or parts of it, is used for research in a
%   publication or if it is fully or partially rewritten for another
%   computing language the authors and institution should be acknowledged
%   in written form and additionally you should cite:
%     M. Liebi, M. Georgiadis, A. Menzel, P. Schneider, J. Kohlbrecher,
%     O. Bunk, and M. Guizar-Sicairos, “Nanostructure surveys of
%     macroscopic specimens by small-angle scattering tensor tomography,”
%     Nature 527, 349-352 (2015).   (doi:10.1038/nature16056)
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its
%    proper use and the correctness of the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%