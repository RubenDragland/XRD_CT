function fourier_figure_fast(proj_show, one_proj, fig)
        processing.load_data = 0; % 1 to force re-load and re-process
        processing.save_fig  = 0;
        processing.plot_data = 1;
        processing.print_fig = 0;
        processing.movie_fig = 0;
        % define which q-range to plot
        processing.which_qrange_to_plot = 1;
        
        four_dat.par=one_proj.par;
        four_dat.integ=one_proj.integ;
        
        four_dat.fnames='';
        
        % only plot the important plots
        four_dat.fmt.fig_sym_int = 0;        % symmetric intensity over all segments: f1_amp
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
        four_dat.fmt.fig_pos_sym_int =  [150 150 screen_size(3)-300 screen_size(4)-300 ];
        four_dat.fmt.fig_pos_comb_all =  [150 150 screen_size(3)-300 screen_size(4)-300 ];
        four_dat.fmt.fig_pos_deg_orient=[];
        four_dat.fmt.fig_pos_orient=[];
        four_dat.fmt.fig_pos_comb_asym=[];
        four_dat.fmt.fig_pos_scatt_asym=[];
        four_dat.fmt.fig_pos_cos_dev=[];
        four_dat.fmt.fig_pos_int_pow_law_exp=[];
        four_dat.fmt.fig_pos_fit_peak_pos=[];
        
        four_dat.fmt.parallel_plotting=false;
        
        four_dat.fmt.cw_size=0.1;
        four_dat.fmt.cw_xoffs=0;
        four_dat.fmt.cw_yoffs=0;
        
        
        four_dat.fmt.interpol=0;
        four_dat.fmt.std_scale=2.0;
        % turn title off
        four_dat.fmt.title_on = 0;
        
        proj_IR = flip(permute(proj_show, [2 3 1]),3);
        
        % calculate the fourier coefficients
        f1_amp = [];
        f2_amp = [];
        f2_phase = [];
        for ii = 1:size(proj_IR, 1)
            
            for jj = 1:size(proj_IR, 2)
                point = squeeze(proj_IR(ii, jj, :));
                [f1_amp1, f2_amp1, f2_phase1] = fourier_coefficients(point);
                f1_amp(ii, jj, :) = f1_amp1;
                f2_amp(ii, jj, :) = f2_amp1;
                f2_phase(ii, jj, :) = f2_phase1;
            end
        end
            
        % Convert the phase of the fft coefficient to orientation in degrees
        % flag for ignoring the sign of f2_phase
        ignore_f2p_sign = 0;% 
%         % plot the simulated symmetric intensity
%         four_dat.fmt.fig_sym_int = 1;        % symmetric intensity over all segments: f1_amp
%         four_dat.fmt.fig_asym_sym_int = 0;   % combination of orientation (colorwheel), degree of orientation (hue)
%         figh = figure(100);
%         subplot(2,3,3)
%         imgh = imagesc(1);
%         four_dat.fmt.figh = figh;
%         four_dat.fmt.imgh = imgh;
%         % plot
%         plot_fourier_result(four_dat, processing);
% %         title(sprintf('IRTT Sym. Int., rotX %.1f, rotY %.1f', one_projection.rot_x, one_projection.rot_y) , 'FontSize', 10);

        phi_det = four_dat.integ.phi_det;
        orientation = phase2orientation(f2_phase,phi_det,ignore_f2p_sign);
        degree_orient = f2_amp ./ (f1_amp + 1e-6); % to avoid a division by zero
        
        % compile return structure for plotting
        four_dat.f1_amp = f1_amp;
        four_dat.f2_amp = f2_amp;
        four_dat.f2_phase = f2_phase;
        four_dat.orientation = orientation;
        four_dat.degree_orient = degree_orient;
% 
%         % plot the simulated symmetric intensity
%         four_dat.fmt.fig_sym_int = 1;        % symmetric intensity over all segments: f1_amp
%         four_dat.fmt.fig_asym_sym_int = 0;   % combination of orientation (colorwheel), degree of orientation (hue)
%         figh = figure(100);
%         subplot(2,3,3)
%         imgh = imagesc(1);
%         four_dat.fmt.figh = figh;
%         four_dat.fmt.imgh = imgh;
%         % plot
%         plot_fourier_result(four_dat, processing);
% %         title(sprintf('IRTT Sym. Int., rotX %.1f, rotY %.1f', one_projection.rot_x, one_projection.rot_y) , 'FontSize', 10);
        
        % plot the simulated orientation
        four_dat.fmt.fig_sym_int = 0;        % symmetric intensity over all segments: f1_amp
        four_dat.fmt.fig_asym_sym_int = 1;   % combination of orientation (colorwheel), degree of orientation (hue)
%         figh = figure(100);
        imgh = imagesc(1);
%         axesh = gca;
        
        four_dat.fmt.figh = fig;
        four_dat.fmt.imgh = imgh;
        % plot
        plot_fourier_result(four_dat, processing);
end