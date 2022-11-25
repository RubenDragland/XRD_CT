%% Generate numerical tensor tomogram phantom
close all;clear;
%addpath(genpath('/home/fredrimu/'));
get(0,'Factory');
set(0,'defaultfigurecolor',[1 1 1])

tomo_sz = [20 20 20 6 ];


a0 = 1; %Define coefficients
a2 = 1/3;
a4 = 1/6;
a6 = 1/10000;
phi_op = pi/2;
theta_op = pi/2;

% Define projections angles
[rot_x,rot_y] = get_ID15_SASTT_angles_new();
%% 


tensorct_phantom = zeros(tomo_sz);

tensorct_phantom(7:13,7:13,7:13,1) = a0; %[Y X Z]

mask3D = squeeze(tensorct_phantom(:,:,:,1)) > 0;
mask3D_lin = mask3D(:);

a0_lin = squeeze(tensorct_phantom(:,:,:,1));
a0_lin = a0_lin(:);

a2_lin = a0_lin.*a2;
a4_lin = a0_lin.*a4;
a6_lin = a0_lin.*a6;

phi_op_lin = ones([1 length(a0_lin)]).*phi_op;
theta_op_lin = ones([1 length(a0_lin)]).*theta_op;
%%
a = [a0_lin' a2_lin' a4_lin' a6_lin'];



%%

load('result_s3t2_reprocessed_09Dec2019_q17-40_all_again_add_id_rcoeff_0.105.mat');



%% generate dummy projection to send into forward proj algorithm
projection = struct();

for ii = 1:length(rot_x)
    projection(ii).data = zeros(30,30,8);
    projection(ii).rot_x = rot_x(ii);
    projection(ii).rot_y = rot_y(ii);
    projection(ii).Rot_exp = calcRotExp(rot_x(ii),rot_y(ii),0,0);
    projection(ii).window_mask = ones(30,30);
       
    
end

params.theta_0 = theta_op_lin;
params.phi_0 = phi_op_lin;
params.a = a;
params.numOfvoxels = length(theta_op_lin);
params.lines = tomo_sz(1);
params.points_per_line = tomo_sz(2);

%%
clear proj_SH;

proj_SH = plot_projection_reg(params,projection(1:259));

%% Build phantom projection from proj_SH

for jj = 1:length(proj_SH)
    
    projection_phantom(jj).data = proj_SH(jj).projection;
    projection_phantom(jj).dx = 0;
    projection_phantom(jj).dy = 0;
    projection_phantom(jj).rot_x = proj_SH(jj).rotx;
    projection_phantom(jj).rot_y = proj_SH(jj).roty;
    projection_phantom(jj).Rot_exp = calcRotExp(proj_SH(jj).rotx,proj_SH(jj).roty,0,0);
    projection_phantom(jj).diode = squeeze(mean(projection_phantom(jj).data,3)).*0.001;
    projection_phantom(jj).par = getParForSASTT('s3t2'); %these are just dummy parameters that follow the data in the cSAXS pipeline
    
    projection_phantom(jj).integ.phi_det = [11.2500   33.7500   56.2500   78.7500  101.2500  123.7500  146.2500  168.7500];
    projection_phantom(jj).integ.norm_sum_avg = ones(8,24);
    projection_phantom(jj).integ.norm_sum = ones(1233,8); 
end

%% Check that data looks okay.

% displayProjections(projection_phantom);

% slighly out of FOV, s3t2
% proj no: 213-216, 233-239, 254-259

% rotates pos y

%% Save data

savename = 'SASTT_test_phantom_7Nov2022_vert_scat.mat';

phantom_params.a0 = a0;
phantom_params.a2 = a2;
phantom_params.a4 = a4;
phantom_params.a6 = a6;
phantom_params.phi_op = phi_op;
phantom_params.theta_op = theta_op;

projection = projection_phantom;

save(savename,'projection','phantom_params');

