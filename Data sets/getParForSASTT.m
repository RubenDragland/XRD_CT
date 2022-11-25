function par = getParForSASTT(sample_name)

       % Copied some dummy paramters from the carbon knot sample for now.

        par.online_analysis = 0;
        par.user_name = 'fkm';
        par.r_sum = {[17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40]};
        par.r_bgr = [0 0;0 0];
        par.I_bgr_pow = [4 4];
        par.peakfit = [0 0 0 1 0];
        par.pix_size= 0.1720;
        par.det_dist = 7126;
        par.lambda_A = 0.9116;
        par.do_segmentation = 0;
        par.pixel_range = 1:2048;
        par.cluster_range = 2;
        par.prepare_SASTT = 1;
        par.which_qrange = 1;
        par.tomo_axis_x = 0;
        par.qresolved_q = [];
        par.ind_max_x = 28;
        par.ind_max_y = 45;
        par.x_scale = 0.0400;
        par.y_scale = 0.0400;
        par.fast_axis_x = 1;
        par.snake_scan = 0;
        
        
        
        % For tensor XRD-CT bone samples.
        
        if strcmp(sample_name,'s8t')
            par.sample_name = 's8t';
            par.ind_max_x   = 65;
            par.ind_max_y   = 65;
        elseif strcmp(sample_name,'s3t2')
            par.sample_name = 's3t2';
            par.ind_max_x   = 65;
            par.ind_max_y   = 67;
        elseif strcmp(sample_name,'shale')
            par.sample_name = 'shale';
            par.ind_max_x = 28;
            par.ind_max_y = 71;
        end
end
            




