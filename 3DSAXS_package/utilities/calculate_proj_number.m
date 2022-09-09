%*------------------------------------------------------------------------*
%|                                                                        |
%|  Except where otherwise noted, this work is licensed under a           |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0             |
%|  International (CC BY-NC-SA 4.0) license.                              |
%|                                                                        |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)     |
%|                                                                        |
%|      Author: CXS group, PSI                                            |
%*------------------------------------------------------------------------*
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sample_diameter = 20; %in micrometers
resolution = 2; %in micrometers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
latitudinal_spacing = 5; %in degrees
number_projections = 360; %only change in case the measurement time is too long for the ideal case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
field_of_view = [30,30]; % in pixels - needed to estimate the duration
exposure_time = 0.5; %in seconds - needed to estimate the duration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N0 = sample_diameter/resolution;
if (isempty(number_projections))
     total_tomo_angle = N0 * 6; % based on the six parameters which should be optimized to solve the spherical harmonics
else
   total_tomo_angle = number_projections;
end
ideal_proj = N0 * 6; % based on the six parameters which should be optimized to solve the spherical harmonics
point_distribution = total_tomo_angle - N0;
tilt_angles = [0:latitudinal_spacing:50];
extra_tilt_angles = length(tilt_angles)-1;
longitudinal_spacing = 180/(point_distribution/extra_tilt_angles);
N_tilt = round((point_distribution/extra_tilt_angles)/2);

number_projections_calc = round(sum((2.*N_tilt.*cosd(tilt_angles))));
overhead = 1.53; %it might change for different stages
readout_time = 0.005; % readout time of the detector in seconds
time = ((field_of_view(1)+1)*(field_of_view(2)+1)*(exposure_time+readout_time)*number_projections_calc*overhead)/3600;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% output: test if that's good enough
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test = latitudinal_spacing/longitudinal_spacing;
if test < 0.8 || test > 1.2
    display(sprintf('\n THE PROJECTIONS ARE NOT WELL DISTRIBUTED. \n Change the latitudinal spacing or the total number of projections.'));
    display(sprintf(' For %d projections (ideal is %d): \n   latitudinal = %.2f degrees \n   longitudinal = %.2f degrees \n   number of extra tilt angles = %d (excluding 0 degrees)\n   use N_tilt = %d \n   estimated duration = %.2f h', number_projections_calc,  ideal_proj, latitudinal_spacing, longitudinal_spacing, extra_tilt_angles, N_tilt, time ));
else
    display(sprintf('\n THE DISTRIBUTION IS OPTIMIZED FOR %d PROJECTIONS (ideal is %d): \n   latitudinal = %.2f degrees \n   longitudinal = %.2f degrees \n   number of extra tilt angles = %d (excluding 0 degrees)\n   use N_tilt = %d \n   estimated duration = %.2f h', number_projections_calc, ideal_proj, latitudinal_spacing, longitudinal_spacing, extra_tilt_angles, N_tilt, time ));
end