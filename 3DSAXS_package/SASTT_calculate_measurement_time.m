%% calculate measurement for SASTT
diameter = 1500;            % diameter of sample in micrometer, this is used calculate number of projections N0
diameter_FOV = 3000;        % diameter of sample for FOV in micrometer, in case something sticks out
height = 2500;%             % height of sample in micrometer
resolution = 45;            % resolution of measurement (equal to beamsize or stepsize, whichever is larger)   
exposure_time = 0.035;      % in sec, including 0.005s for readout
overhead_per_line = 3.4;    % overhead per line, typically 3-4s
overhead_rot = 8;           % overhead per projection
overhead_tilt = 75;         % overhead per tilt and radiation damage check

%%% Change the following parameters to adjust measurement time %%%
num_tilts_reduced = 5;      % number of tilt angles, can be reduced for time reasons
N0_tilt = 20;               % number of projections per tilt angle (not for tilt = 0). 
                            % Ideal sampling when N0_tilt = N0, but usually can be reduced for time reasons

%%
time_per_proj = (exposure_time * (height/resolution +1)+ overhead_per_line) * (diameter_FOV/resolution+1);
time_per_proj_min = time_per_proj/60;

N0 = floor(diameter/resolution);   % number of projections at tilt 0
Ntot1 = 6.*N0;                      % based on reconstructing 6 parameters per voxel
deltaalpha = 180/N0;              % stepsize of rotation motor deltaalpha
deltaalpha_reduced = 180/N0_tilt;
max_tilt = 45;                      % maximum tilt angle available, typically 45deg
num_tilts_ideal = max_tilt/deltaalpha; % nr of tilts having same angular sampling so deltaalpha = deltabeta


deltabeta_ideal = max_tilt/num_tilts_ideal;
tilt_angles_ideal = [0:deltabeta_ideal:max_tilt];
num_projections_tilt_ideal = 2*N0.*cosd(tilt_angles_ideal);
angularstep_intilt_ideal = 360./num_projections_tilt_ideal;
Ntot_ideal = N0 + sum(num_projections_tilt_ideal(2:end));


deltabeta_reduced = max_tilt/num_tilts_reduced;
tilt_angles_reduced = [0:deltabeta_reduced:max_tilt];
num_projections_tilt_reduced = 2*N0_tilt.*cosd(tilt_angles_reduced);
angularstep_intilt_reduced = 360./num_projections_tilt_reduced;
Ntot_reduced = N0 + sum(num_projections_tilt_reduced(2:end));

%%
measurement_time = ((time_per_proj+overhead_rot) * floor(Ntot_reduced) + (overhead_tilt+time_per_proj)*num_tilts_reduced)/3600;

measurement_time_ideal_sampling = (time_per_proj * floor(Ntot_ideal))/3600;

fprintf('********************************\n');
fprintf('resolution: %d um diameter: %0.1f mm height: %0.1f mm h\n',resolution,diameter_FOV/1000,height/1000)
fprintf('measurement time         %0.1f h     ideal sampling %0.1f h\n', measurement_time, measurement_time_ideal_sampling)
fprintf('Number of projections    %d        ideal sampling %d        Based on reconstructing 6 parameters per voxel:  %d\n', floor(Ntot_reduced),floor(Ntot_ideal),floor(Ntot1))
fprintf('delta_alpha at tilt 0    %0.1f deg    ideal sampling %0.1f deg\n', deltaalpha, deltaalpha)
fprintf('delta_alpha              %0.1f deg    ideal sampling %0.1f deg\n', deltaalpha_reduced, deltaalpha)
fprintf('delta_beta               %0.1f deg    ideal sampling %0.1f deg\n', deltabeta_reduced, deltabeta_ideal)
fprintf('*** Input for SASTT template ***\n');
fprintf('N0 = %d   N0_tilt = %d\n', N0,N0_tilt)
fprintf('tilt angles %0.1f\n', tilt_angles_reduced)
