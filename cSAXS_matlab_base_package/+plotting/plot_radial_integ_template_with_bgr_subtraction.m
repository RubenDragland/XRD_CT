% plot_radial_integ_template_with_bgr_subtraction.m
%
% This is a template like plot_radial_integ_template.m which does
% background subtraction by using the plot_radial_integ.m to do averaging
% of segments and/or points and returning the q, intensity vectors, which
% are then plotted afterwards.
% 
% With this template transmission corrections can be done from SPEC
% counters like diode or cyb (mcs not yet supported)

scan_nr=2795;         
backgroundscan=2957 ;  % leave empty [] for no background correction. If a scan is specified
                      % Careful: if not empty, make sure that parameters in
                      % plot_radial_integ bellow are the same for both
                      % measurement and background
point_range=[];       % Here points refer to points in a scan excluding burst frames. For cont_line scans point_range=1; 
                      % Leave empty for all points. Example: point_range=[1:5]
burst_range=[56];       % Subset of burst frames. If empty [] it takes all burst frames                       

base_path='~/Data10/';
eaccount=beamline.identify_eaccount;
detno = 1;                  % Detector number (1 for Pilatus 2M, 2 for Pilatus 300k, 3 for Eiger).
                            % If different than 1, file_flag should be changed e.g. to '_waxs' or '_waxs_2M_at_two_meters'
file_flag='';               % Example: file_flag='_2M_at_two_meters' (or '_waxs' or '_waxs_2M_at_two_meters'  for detno = 2) 
                            % Leave empty '' for default names
transmission_det='diode';   % Leave empty '' for no tranmission correction 
                            % Currently only spec value supported (no mcs)
                            % Only needed for background subtraction

args = {...
'NewFig',0,...              open a new figure for each file, default is 0
'ClearFig',1,...            clear the figure before plotting, default is 1
'XLog',1,...                logarithmic scaling of the x-axis, default is 0
'YLog',1,...                logarithmic scaling of the y-axis, default is 1
'PlotQ',0,...               plot as a function of momentum transfer q rather than pixel no., default is 0
'PlotAngle',0,...           plot as a function of the azimuthal angle rather than q or the radius, default is 0
'RadiusRange',[],...        for azimuthal plots the intensity over this radius range is averaged, default is [] for all radii
'QMulPow',[],...            multiply intensity with q to the power of this value, default is [ ] for no multiplication
'Inverse_nm',1,...          plot q in inverse nm rather than inverse Angstroem, default is 0
'SegAvg',1,...              average over angular segments rather than plotting them with different line colours, default is 1
'SegRange',[],...           segment range to plot, default is [] for all segments
'LegendMulSeg',1,...        show a legend in case of multiple segments being plotted, default is 1
'PointAvg',1 ...%            %plot the average of all intensity curves in the file, which typically means the average of a scan line, default is 1
};

spec_data = io.spec_read(base_path,'ScanNr',scan_nr);
burst_frames=spec_data.burstn;

if isempty(burst_range)
    subset_frames = 1:burst_frames;
else
   if numel(burst_range) > burst_frames
       error(sprintf('The number of elements of burst_range is %d. It must be equal or smaller than burst_at_each_point=%d ',numel(burst_range),burst_frames))
   elseif max(burst_range) > burst_frames
       error(sprintf('The maximum value of burst_range is %d. It must be equal or smaller than burst_at_each_point=%d',max(burst_range),burst_frames))
   else    
       subset_frames=burst_range;
   end
end

if burst_frames ~= 1
    point_range_frames=zeros(1,numel(subset_frames)*numel(point_range));
    hh=1;
    for ii=point_range
        for jj=subset_frames
            point_range_frames(hh)=burst_frames*(ii-1)+jj;
            hh=hh+1;
        end   
    end
elseif burst_frames == 1
    point_range_frames=point_range;
end

if isempty(backgroundscan) == 0
    spec_data_bgr = io.spec_read(base_path,'ScanNr',backgroundscan);
    if isempty(transmission_det)
        correction = 1;
    else
        trans_data=spec_data.(transmission_det);
        trans_data_bgr=spec_data_bgr.(transmission_det);        
        if isempty(point_range)
            transmission_data=mean(trans_data);
            transmission_bgr=mean(trans_data_bgr);
        else
            transmission_data=mean(trans_data(point_range));
            transmission_bgr=mean(trans_data_bgr(point_range));
        end
        correction=transmission_data/transmission_bgr;
    end
end

args={args{:},'PointRange',point_range_frames};

datafile=fullfile(base_path,'analysis',sprintf('radial_integration_%s', file_flag),sprintf('%s_%d_%05d_00000_00000_integ.mat',eaccount,detno,scan_nr));
backgroundfile=fullfile(base_path,'analysis',sprintf('radial_integration_%s',file_flag),sprintf('%s_%d_%05d_00000_00000_integ.mat',eaccount,detno,backgroundscan));

args1={args{:},'FigNo',1};
[x,y]=plotting.plot_radial_integ(datafile, args1);

if isempty(backgroundscan) == 0

args2={args{:},'FigNo',2};
[xb,yb]=plotting.plot_radial_integ(backgroundfile, args2);

figure(3);
loglog(x,y); hold on
loglog(xb,correction*yb); hold off
legend(sprintf('Scan S%05d',scan_nr),sprintf('Scan S%05d',backgroundscan))

figure(4);
loglog(x,y-correction*yb);
title(sprintf('Scan S%05d - S%05d',scan_nr,backgroundscan))

end

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
