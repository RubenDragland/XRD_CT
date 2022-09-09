% Template to plot multi-channel scaler (MCS) data for scanning based on
% the frameID of scanning SAXS

base_path = '~/Data10';
channel_to_plot = 3;

% To plot all the frames starting from [starting_frame], put
% plot_all_frames=1;
plot_all_frames = 0;
starting_frame  = 2;
wait_for_next   = 0; % If to wait until next frame starts

% Or choose the frame_ids you want to plot (used if plot_all_frames = false):
Frames = [1:10];


%% Modify until here

if (plot_all_frames)
    i=starting_frame;
    while 1==1
        [scann,p] = get_scannum_frame(base_path,i);
        while (isempty(scann))
            fprintf('Frame %d not started, waiting 10 secs...\n',i);
            pause(10);
            [scann,p] = get_scannum_frame(base_path,i);
        end
        if (wait_for_next == 1)
            [scann,p] = get_scannum_frame(base_path,i+1);
            while (isempty(scann))
                fprintf('Next Frame not started, waiting 10 secs...\n');
                pause(10);
                [scann,p] = get_scannum_frame(base_path,i+1);
            end
        end
        try
            [scann,p] = get_scannum_frame(base_path,i);
            if (p.fast_axis_y == 1)
                nscan = p.saxs_intervals(1);
            else
                nscan = p.saxs_intervals(2);
            end
            io.mcs_mesh(scann(1), nscan ,'ChToPlot',channel_to_plot, 'FastAxisX', 1-p.fast_axis_y , 'SnakeScan', p.snake_scan, 'Pos_file', base_path, 'FrameID', p.frame_id);
        catch err
            fprintf('Error occurred with frame %d, skipping it...\n',i);
        end
        i=i+1;
    end
else    
    for i=1:numel(Frames)
        [scann,p]=get_scannum_frame(base_path,Frames(i));
        io.mcs_mesh(scann(1), scann(end)-scann(1) ,'ChToPlot',channel_to_plot, 'FastAxisX', 1-p.fast_axis_y , 'SnakeScan', p.snake_scan, 'Pos_file', base_path, 'FrameID', p.frame_id);
    end
end

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
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
%   using the “cSAXS matlab package” developed by the CXS group,
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
