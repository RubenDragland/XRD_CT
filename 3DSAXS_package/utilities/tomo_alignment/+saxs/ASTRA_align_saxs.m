%   FUNCTION [optimal_shift,par, rec, range] = ...
%             ASTRA_align_tomo(sinogram_0,weights, angles, Npix, optimal_shift, par, varargin )
%       self consitency based alignment procedure based on the ASTRA toolbox 
% Inputs: 
%       sinogram_0 - real value sinogram (ie not diff)
%       angles - angle in degress
%       Npix - size of the reconstructed field 
%       optimal_shift - initial guess of the shift 
%       par - parameter structure 
%           -> params:
%                   align_vertical - allow vertical alignement 
%                   align_horizontal - allow horizontal alignement 
%                   mask_threshold - values < threshold will be considered
%                       empty, [] == auto guess 
%                   apply_positivity - remove negatie values 
%                   high_pass_filter - bandpath filter applied on the data
%                   to avoid low spatial freq. errors 
%                   max_iter - maximal number of iterations for the
%                   alignment 
%                   lamino_angle - tilt angle of the tomographic axis 
%                   show_diff_movie - At the end of alignment it shows a
%                       movie with the projections and the difference
%                   valid_angles - Binary array for determining which projections are used in the reconstruction of the 
%                       tomogram used for alignment. The non-valid projections are not used for the tomogram but they are aligned, 
%                       useful when 45 degrees gives artifacts. 



function [optimal_shift,par, rec, sinogram_shifted] = ...
            ASTRA_align_saxs(sinogram_0, angles, Npix, optimal_shift, par, varargin )
    


   % standard values: usually not needed to change for users
    par.center_reconstruction = true; % keep the center of mass in center of rec. volume
    par.lamino_angle = 90 + par.tilt_angle;  % deg - laminography angle
    Nprojections = length(angles);
    par.use_mask = 0;    % apply support mask
    par.mask_threshold = []; % empty == Otsu thresholding
    par.apply_positivity = 1; % remove negative values

    parser = inputParser;
    parser.addParameter('align_vertical',  true , @islogical )
    parser.addParameter('align_horizontal', 0.5 , @isnumeric )
    parser.addParameter('mask_threshold',  [] , @isnumeric )
    parser.addParameter('use_mask',  false , @islogical )
    parser.addParameter('apply_positivity',  true , @islogical )
    parser.addParameter('high_pass_filter',  0.02 , @isnumeric )
    parser.addParameter('max_iter',  50 , @isint )
    parser.addParameter('lamino_angle',  90 , @isnumeric )
    parser.addParameter('show_diff_movie',  0 , @islogical )
    parser.addParameter('valid_angles',  1:Nprojections > 0 , @islogical )

    parser.parse(varargin{:})
    r = parser.Results;

    % load all to the param structure 
    for name = fieldnames(r)'
        if ~isfield(par, name{1})  % prefer values in param structure 
            par.(name{1}) = r.(name{1});
        end
    end

    
    import saxs.*
        
    [Nlayers,~,Nangles] = size(sinogram_0);

if par.smooth_data
   sinogram_0 =imgaussfilt(sinogram_0, par.smooth_data); 
end
    
%     disp('Shifting sinograms')
    sinogram_0 = smooth_edges(sinogram_0); 
    % shift to the last optimal position 
%     for i = 1:Nangles
%         sinogram_0(:,:,i) =  imshift_fast(sinogram_0(:,:,i), -optimal_shift(i,1), -optimal_shift(i,2), [], 'nearest');
%     end
    sinogram_0 =  utils.imshift_fft(sinogram_0, optimal_shift);
    
    %% limit the vertical range 
%     disp('Cropping vertical range')
    range =  [max(par.vert_range(1),round(2+max(optimal_shift(:,2)))),min(par.vert_range(end), floor(Nlayers-1 + min(optimal_shift(:,2))))]; 
    if (isempty(range)||range(1)>=range(2)) 
        error('Too small par.vert_range'); 
    end 
    sinogram_0 = sinogram_0(range(1):range(2),:,:); 
    sinogram_0 = smooth_edges(sinogram_0); 
    
        
                
    %% %%%%%%%%%%%%%%%% initialize astra %%%%%%%%%%%%%%%%
    [Nlayers,width_sinogram,~] = size(sinogram_0);
      
    
    min_value= 0;
    MAX_SINO = math.sp_quantile(sinogram_0, 0.999,10); 
    MASS = math.mean2(sinogram_0); 
    W = tukeywin(width_sinogram, 0.2)';  % avoid edge issues 
    if Nlayers > 10 && par.align_vertical ; W = tukeywin(Nlayers, 0.2)*W; end 
    clear err 
    
    lamino_angle = par.lamino_angle(:); 
    
    shift_all = nan(par.max_iter,Nangles,2);
    shift_all(1,:,:) = 0;

    time_plot = tic; 
    
    
    for ii = 1:par.max_iter-1
        fprintf('Iteration %i \n ', ii)
        %% shift the sinogram
        sinogram_shifted =  utils.imshift_fft(sinogram_0, squeeze(shift_all(ii,:,1)), squeeze(shift_all(ii,:,2)));

        % get configuration for ASTRA, may differ in every iteration 
        [cfg, vectors] = ...
            ASTRA_initialize(Npix, Nlayers,width_sinogram,angles,lamino_angle,0, 1); 
        % find optimal split of the dataset for given GPU 
        split = 1;

        %% FBP method 
        rec  = FBP(sinogram_shifted, cfg, vectors, split,par.valid_angles);
       
%         keyboard
        
        
        
        %% %%%%%%%%% apply all availible prior knowledge %%%%%%%%%%%%%%%%%%%%%%%
        if par.use_mask
            if mod(ii,5) == 1 
                %% get mask 
                if isempty(par.mask_threshold)
                    mask = rec > graythresh(rec(:));
                else
                    mask = rec > par.mask_threshold; 
                end
                disp('Find mask')
                for i = 1:Nlayers
                    utils.progressbar(i, Nlayers)
                    mask(:,:,i) = imfill(imdilate(mask(:,:,i), strel('disk', 1)),'holes');
                end
                figure(46)
                subplot(1,2,1)
                imagesc(rec(:,:,ceil(end/2)) .* ~mask(:,:,ceil(end/2)))
                colorbar
                axis off image 
                colormap bone 
                title('Example of residuum after applied mask')
                subplot(1,2,2)
                hist(rec(1:100:end), 100)
                axis tight 
                drawnow 
                
                min_value = median(rec(1:100:end));
                               
                
            end
            %% apply mask 
            rec = rec .* mask;
        end
        
        if par.apply_positivity
            rec = max(0,rec);
            if min_value > 0
                % assume minimal value inside mask regions 
                rec = rec.*~mask + mask.*max(min_value,rec);
            end
        end
                

        if par.center_reconstruction && ii < 20
            %% try to keep reconstruction in center 
            [x,y,mass] = math.center(sqrt(max(0,rec))); 
            % more robust estimation of center 
            indx = find(~isnan(squeeze(x)));
            indy = find(~isnan(squeeze(y)));
            x = mean(x(indx).*mass(indx).^2)./mean(mass(indx).^2); 
            y = mean(y(indy).*mass(indy).^2)./mean(mass(indy).^2); 
            rec = utils.imshift_fft(rec, -x/5, -y/5); % go slowly
        end

        %% %%%%%%%%  Get computed sinogram 
        disp('Forward projection')
        % forward projection 

        sinogram_corr = saxs.Ax_CPU(rec,cfg, vectors);
         
        %% find shift of the data sinogram to match the "optimal sinogram"
        disp('Find optimal shift')
       
        %% find optimal shift 
        [dX,dY] = get_img_grad(smooth_edges(sinogram_corr), 1);
        DS = sinogram_corr - sinogram_shifted;
                
        %% apply high pass filter => get rid of phase artefacts 
        DS = imfilter_high_pass_1d(DS,2,par.high_pass_filter); 
        if par.align_horizontal
            % calculate optimal shift of the 2D projections in horiz direction 
            dX = imfilter_high_pass_1d(dX,2,par.high_pass_filter); 
            shift_X = -squeeze(math.sum2(W.* dX .* DS) ./ math.sum2(W.*dX.^2));
        end
        if par.align_vertical    
            % calculate optimal shift of the 2D projections in vertical direction 
            dY = imfilter_high_pass_1d(dY,2,par.high_pass_filter); 
            shift_Y = -squeeze(math.sum2(W.* dY .* DS) ./ math.sum2(W.* dY.^2));
            shift_Y = shift_Y - mean(shift_Y);
        else
            shift_Y = zeros(Nangles,1); 
        end
            
%         keyboard
        
        
        err(ii,:) = gather(sqrt(math.mean2(DS.^2)) ./ MASS); 
        
        if any(isnan(shift_X)) || any(isnan(shift_Y))
            warning('Shift is nan')
            keyboard
        end
               

        % do not allow more than 1px per iteration !! 
        shift_X = min(1, abs(shift_X)).*sign(shift_X);
        shift_Y = min(1, abs(shift_Y)).*sign(shift_Y);
                
        step_relaxation = 0.5; %% small relaxation is needed to avoid oscilations  
        shift_all(ii+1,:,:) = shift_all(ii,:,:) + gather(reshape([shift_X, shift_Y], [1,Nangles,2]))*step_relaxation;
        % remove degree of freedom in the vertical dimension (avoid drifts)
        shift_all(ii+1,:,2) = shift_all(ii+1,:,2) - median(shift_all(ii+1,:,2));
        
        fprintf('Maximal step update: %3.2g px stopping criterion: %3.2g px\n \n ', max(max(abs(shift_X(par.valid_angles))),max(abs(shift_Y(par.valid_angles)))), par.min_step_size)
        if max(max(abs(shift_X(par.valid_angles))),max(abs(shift_Y(par.valid_angles)))) < par.min_step_size
            break  % stop iterating if converged 
        end
        
        if toc(time_plot) > 1  % avoid plotting in every iteration, it is too slow ... 
        %%%%%%%%%%%%%%%%%%%% plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        disp('plotting')
        figure(5464)
        clf()
        subplot(2,3,1)
        imagesc(squeeze(sinogram_0(ceil(Nlayers/2),:,:))', [0,MAX_SINO]);  
        axis off 
        colormap jet 
        title('Original sinogram')

        subplot(2,3,2)
        hold on 
        plot(shift_X, 'r')
        plot(shift_Y, 'b')
        plot(find(par.valid_angles), err(ii,par.valid_angles)*10, 'k.') 
        plot(find(~par.valid_angles),err(ii,~par.valid_angles)*10, 'r.') 
        hold off 
        legend({'horiz', 'vert', 'errors'}, 'Location', 'northwest')
        title('Current position correction')
        xlim([1, Nangles])
        grid on
        subplot(2,3,3)
        imagesc(squeeze(sinogram_shifted(ceil(Nlayers/2),:,:))', [0,MAX_SINO]);  
        axis off 
        colormap bone 
        title('Corrected sinogram')
        subplot(2,3,4)
        try
        slice = rec(:,:,ceil(Nlayers/2));
        slice = rot90(slice);
        imagesc(slice, [0, math.sp_quantile(rec(:,:,ceil(Nlayers/2)), 0.999,1)])
        xlabel('x')
        ylabel('z')
        axis equal tight xy
        end
        %axis off image
        title(sprintf('Reconstruction, slice=%.0d', ceil(Nlayers/2)+par.vert_range(1)  ))
        colormap bone 
        subplot(2,3,5)
        hold on
        plot(err)
        plot(mean(err,2), 'k', 'LineWidth', 3);
        hold off 
         grid on 
         axis tight
        xlim([1,ii+1])
        set(gca, 'xscale', 'log')
        set(gca, 'yscale', 'log')
        title('Tomogram mispositioning')
        xlabel('Iteration')
        ylabel('Mean abs error')
        
        subplot(2,3,6)
        tmp = [shift_all(:,:,1), ones(par.max_iter,5)*max(err(:)), shift_all(:,:,2)]; 
        imagesc(tmp)
        colorbar 
        ylabel('Iteration')
        xlabel('Projection id')
        colormap bone 
        title('dx shift        dy shift')
%         print('-dpng', sprintf('FBP_pos_grad_solver_3D_%02i.png', ii))
        drawnow   
        
        time_plot = tic; 

    end
    end
            
    %% prepare outputs to be exported     
    % sort / crop the angles 
    optimal_shift = optimal_shift + squeeze(shift_all(ii+1,:,:));

    % center only in the first  binning round 
    par.center_reconstruction = false; 
    disp('Done')
    %% Show movie of alignment
    if par.show_diff_movie
        figure(1000)
        clf
        subplot(2,3,1);
        colormap bone
        axis xy equal tight
        h1 = imagesc(sinogram_shifted(:,:,1));
        
        subplot(2,3,2);
        colormap bone
        axis xy equal tight
        colorbar
        h2 = imagesc(sinogram_corr(:,:,1));
        
        
        subplot(2,3,3);
        colormap bone
        axis xy equal tight
        colorbar
        rms_error = sqrt(((sinogram_corr-sinogram_shifted).^2)./numel(sinogram_shifted(:,:,1)));
        h3 = imagesc(rms_error(:,:,1));
       
        
        
        if isempty(par.which_projection)
            index = 1:par.stepping:size(sinogram_shifted,3);
        else
            index = par.which_projection;
        end
        
        for ii = index
            h1.CData = sinogram_shifted(:,:,ii);
            h2.CData = sinogram_corr(:,:,ii);
            h3.CData = rms_error(:,:,ii);
            subplot(2,3,1); axis xy equal tight; colorbar
            title(sprintf('Shifted projection \n proj. %d', ii))
            subplot(2,3,2); axis xy equal tight; colorbar
            title(sprintf('Computed projection\n proj. %d', ii))
            subplot(2,3,3); axis xy equal tight; colorbar
            title(sprintf('RMS error\n proj. %d', ii))
            
            subplot(2,3,[4:6]);
            errorbar(ii, mean(mean(rms_error(:,:,ii).^2)),std(std(rms_error(:,:,ii).^2)), 'ok')
            hold on
            grid on
            xlabel('Projection')
            ylabel('mean squared error')
            axis tight
            drawnow
        end
    end
    
end