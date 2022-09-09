% Adjust projections - With this script you can do some simple operations
% on the projections. For example fixing offsets between lines or cropping
% the projections

addpath ../
% import utils.*
clear
flag_adjust_positions = true;
flag_resample_y = false;  % if = 3 the y axis will have 3 times more pixels
proj_filename   = '/das/work/units/csaxs/p17917/analysis/SASTT/cortical_S2/projection_data/SASTT_cortical_S2.mat';
extrasuffix     = '_corrected_2';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
proj     = load(proj_filename);
proj_new = proj;
[thepath,thename,theext] = fileparts(proj_filename);
proj_filename_new = fullfile(thepath,[thename extrasuffix theext]);



proj_fields = fieldnames(proj);

data_plot_new = zeros([size(proj.(proj_fields{1}).data,1) size(proj.(proj_fields{1}).data,2) numel(proj_fields)]);
data_plot_old = data_plot_new;


for ii = 1:numel(proj_fields)
    fprintf('Adjusting SASTT projection %s\n', proj_fields{ii})
    current_proj = proj.(proj_fields{ii}); 
    current_proj_before = current_proj;
    
    
    if flag_adjust_positions
        avy = diff(mean(proj.(proj_fields{ii}).positions_out(:,:,2),1));
        avy = -mean(avy.*(-1).^([1:numel(avy)]));
        avy_step = -mean(mean(diff(proj.(proj_fields{ii}).positions_out(:,:,2),[],1)));
        snake_scan_pos_corr = avy/avy_step; % Correction of each second line
        
        fprintf('Correcting snake scan position offset in %s \n',proj_fields{ii})
        fprintf('    Average y step size %d microns \n',avy_step*1e3)
        fprintf('    Average y error     %d microns \n',avy*1e3)
        fprintf('    Correction of       %d pixels \n',snake_scan_pos_corr)
        fprintf('    --------------------------- \n')
        
        
        for jj = 1:size(current_proj.data,3)
            current_proj.data(:,2:2:end,jj) = utils.shiftwrapbilinear(current_proj.data(:,2:2:end,jj),snake_scan_pos_corr,0);
        end
        current_proj.diode(:,2:2:end)           = utils.shiftwrapbilinear(current_proj.diode(:,2:2:end),snake_scan_pos_corr,0);
        current_proj.scan_num(:,2:2:end)        = utils.shiftwrapbilinear(current_proj.scan_num(:,2:2:end),snake_scan_pos_corr,0);
        current_proj.scan_point(:,2:2:end)      = utils.shiftwrapbilinear(current_proj.scan_point(:,2:2:end),snake_scan_pos_corr,0);
        current_proj.positions_out(:,2:2:end,2)   = current_proj.positions_out(:,2:2:end,2) - avy;
        
    end
    
    if flag_resample_y
%         outputsizey = flag_resample_y*size(current_proj.diode,1);
%         aux = real(utils.interpolateFT_ax(current_proj.diode,outputsizey,1));
        
        outputsize = size(current_proj.diode);
        outputsize(1) = outputsize(1)*3;
        current_proj.par.ind_max_y    = outputsize(1);
        current_proj.par.y_scale      = proj_new.(proj_fields{ii}).par.y_scale/flag_resample_y;
        
        outputsizedata = size(current_proj.data);
        outputsizedata(1) = outputsizedata(1)*3;

        current_proj.diode      = imresize(current_proj.diode,outputsize,'bilinear');
        if isfield(current_proj,'scan_num')
            current_proj.scan_num   = imresize(current_proj.scan_num,outputsize,'bilinear');
        end
        if isfield(current_proj,'scan_point')
            current_proj.scan_point = imresize(current_proj.scan_point,outputsize,'bilinear');
        end
        
        current_proj.positions_out = [];
        current_proj.data = imresize3(current_proj.data, outputsizedata, 'linear');
        
    end
    
%     warning('Make this more general, mising diode, scan numbers, positions, etc')
    proj_new.(proj_fields{ii}) = current_proj;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Prepare for plotting %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    data_plot_new(:,:,ii) = current_proj.data(:,:,1);
    data_plot_old(:,:,ii) = current_proj_before.data(:,:,1);
end



%%
figure(1);
subplot(1,2,1)
title('Before correction')
plotting.imagesc3D(data_plot_old); 
axis xy equal tight
colormap bone
drawnow
colorbar
subplot(1,2,2)
title('After correction')
plotting.imagesc3D(data_plot_new); 
axis xy equal tight
colormap bone
drawnow
colorbar

[filepath, filename, fileext] = fileparts(proj_filename);
filename = fullfile(filepath,[filename '_adjusted' fileext]);

    
save(proj_filename_new, '-struct', 'proj_new')
display(['Saving SASTT projection to: ' proj_filename_new])


%% Plot one set of positions
position_choose = 52;
x = proj_new.(proj_fields{position_choose}).positions_out(:,:,1);
y = proj_new.(proj_fields{position_choose}).positions_out(:,:,2);
figure(2); 
subplot(1,2,2)
plot(x(:),y(:),'o');
title('New positions')
axis equal

x_old = proj.(proj_fields{position_choose}).positions_out(:,:,1);
y_old = proj.(proj_fields{position_choose}).positions_out(:,:,2);
subplot(1,2,1)
plot(x_old(:),y_old(:),'o');
title('Old positions')
axis equal


