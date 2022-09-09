%ONLINE_VIEWER plot last scan 
% 
% ** basePath           path to the data directory; default: ~/Data10
% ** dataDir            name of the data directory; default: data
% ** filter             filter for nexus_read
%
%
% EXAMPLE:
%       plotting.online_viewer();
%       plotting.online_viewer('filter', 'pilatus_2');
%       plotting.online_viewer('filter', 'eiger_4', 'mask', '~/Data10/cxs_software/base/detector_masks/eiger_4_valid_mask.h5');
%       
%
% see also: io.nexus_read()

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


function online_viewer(varargin)

par = inputParser;
par.addParameter('dataDir', 'data', @ischar)
par.addParameter('filter', '/entry/data', @ischar)
par.addParameter('basePath', '~/Data10', @ischar)
par.addParameter('useLog', true, @islogical)
par.addParameter('addProfiles', false, @islogical)
par.addParameter('mask', [], @ischar)

par.parse(varargin{:})
vars = par.Results;

prev_scanNr = -1;
scanNr = -1;




while true
    try
        % get parent directory (S00000-00999)
        parentDir = dir(fullfile(utils.abspath(vars.basePath), vars.dataDir));
        if length(parentDir(end).name)>3
            % get scan directory (S00250)
            subDir = dir(fullfile(parentDir(end).folder,  parentDir(end).name));
            if length(subDir(end).name)>3
                scanNr = str2double(subDir(end).name(2:end));
                if scanNr ~= prev_scanNr
                    % get file name
                    fname = utils.find_nexus_file(fullfile(utils.abspath(vars.basePath), vars.dataDir), scanNr);
                    waitTime = 0;
                    fready = true;
                    while ~utils.nexus_file_is_ready(fname)
                        pause(1)
                        waitTime = waitTime + 1;
                        if waitTime > 20
                            fready = false;
                            break;
                        end
                    end
                    
                    if fready
                        % read data
                        args = {};
                        if ~isempty(vars.filter)
                            args = [args{:} {'filter', vars.filter}];
                        end
                        if ~isempty(vars.mask)
                            args = [args{:} {'mask', vars.mask}];
                        end
                        [data, m] = io.nexus_read(fname, args{:});

                        if isempty(data)
                            warning('Did not find %s in %s',vars.filter,fname)
                            pause(1)
                            continue
                        end
                        if ~isempty(vars.mask)

                            data.data = data.data.*repmat(m,[1 1 size(data.data,3)]);
                        end
                        if vars.useLog
                            data.data = log10(abs(single(data.data)+1e-10));
                        end

                        int1 = squeeze(sum(data.data,1));
                        int2 = squeeze(sum(data.data,2));
                        
                        % plot
                        color_axis = math.sp_quantile(data.data, [1e-6 1-1e-6], 10);
                        color_axis(1) = max([color_axis(1) 0]);
                        if vars.addProfiles
                            ax = {};
                            fig = plotting.smart_figure(1258);
                            clf();
                            subplot(6,6,[1 5]);
                            ax{end+1} = gca();
                            subplot(6,6,[[1:6]*6]);
                            ax{end+1} = gca();
                            vals = reshape([1:36], 6, 6)';
                            vals(1,:) = [];
                            vals(:,end) = [];
                            subplot(6, 6, [vals(:)]);
                            ax_controller = gca();
                            
                            
                            plotting.imagesc3D(data.data); axis equal tight xy;
                            caxis(color_axis)
                            colormap(plotting.colormaps.franzmap);
                            colorbar
                            plot(ax{1}, int1(:,1)); axis equal tight;
                            set(ax{1},'xtick',[])
                            set(ax{1},'xticklabel',[])
                            xlim(ax{1}, [1 size(int1(:,1),1)]);
                            ax{1}.Position = [ax_controller.Position(1) ax_controller.Position(2) + ax_controller.Position(4)+0.01 ax_controller.Position(3)  0.15];
                            plot(ax{2}, int2(:,1)); view(ax{2}, [90 -90]); axis equal tight
                            set(ax{2},'xtick',[])
                            set(ax{2},'xticklabel',[])
                            xlim(ax{2}, [1 size(int2(:,1),1)]);
                            ax{2}.Position = [ax_controller.Position(1)+ax_controller.Position(3)+0.01 ax_controller.Position(2) 0.15 ax_controller.Position(4)];
                            if ~isempty(vars.filter)
                                plotting.suptitle(sprintf('S%05d - %s', scanNr, vars.filter))
                            else
                                plotting.suptitle(sprintf('S%05d', scanNr))
                            end
                            if size(data.data,3)>1
                                propListener = addlistener(ax_controller.slider_handle,'Value','PostSet',@(src,evnt)change_profiles(ax_controller, int1, int2, ax));
                                ax_controller.play(ax_controller);
                            end
                            
                        else
                            fig = plotting.smart_figure(1258);
                            plotting.imagesc3D(data.data); axis equal tight xy;
                            caxis(color_axis)
                            colormap(plotting.colormaps.franzmap);
                            colorbar
                            ax_controller = gca();
                            if ~isempty(vars.filter)
                                title(sprintf('S%05d - %s', scanNr, vars.filter))
                            else
                                title(sprintf('S%05d', scanNr))
                            end
                            
                            if size(data.data,3)>1
                                ax_controller.play(ax_controller);
                            end
                        end


                        prev_scanNr = scanNr;
                    end
                    
                end
                pause(1)
            end
        end
    catch ME
        keyboard
    end
    
end




end



function change_profiles(ax_controller, int1, int2, ax)
set(ax{1}.Children, 'ydata', int1(:,round(ax_controller.slider_handle.Value)))
set(ax{2}.Children, 'ydata', int2(:,round(ax_controller.slider_handle.Value)))
ax{2}.Position = [ax_controller.Position(1)+ax_controller.Position(3)+0.01 ax_controller.Position(2) 0.15 ax_controller.Position(4)];
ax{1}.Position = [ax_controller.Position(1) ax_controller.Position(2) + ax_controller.Position(4)+0.01 ax_controller.Position(3)  ax{2}.Position(3)];

end


