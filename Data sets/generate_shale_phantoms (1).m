%% Generates Shale phantom
close all;clear;
addpath(genpath('G:/PhD/Work/codebase/'));
get(0,'Factory');
set(0,'defaultfigurecolor',[1 1 1])

%%
% The shale phantom is a cylinder/disk with orientation pointning mainly
% in one direction but with a distribution
% This script generates N phantoms where the main local voxel orientation
% is tilted with the y (tomo) axis, ranging from 0 degree to 90 degrees.

%% Load the original shale data from clin001
% load([pwd '/SASTT_shale_q_38.mat']); %loads projections
load([pwd '/SASTT_PS1_peak_no_1_bgsub.mat']);


% load([pwd '/result_shale_q_38_q17-40_all_again_test_cap_reg_coeff_0.000.mat']);
load([pwd '/result_PS1_clin001_q17-40_all_again_filloutliers_reg_coeff_5.000.mat']);
%loads reconstructed data

%% Compare measured shale data to a simulated phantom to get a guess of the a's.
% f_prev_sim_ph = load([pwd '/shale_phantom_28Aug2020_phi_op_mean_1.57_aligned_ASTRA.mat']); %loads previously simulated phantom.
% f_real_data = load([pwd '/SASTT_shale_q_38.mat']);
% point = [40 19];
% figure(111);
% clf
% hold on;
% plot(squeeze(f_real_data.projection(1).data(point(2),point(1),:)));
% plot(squeeze(f_prev_sim_ph.projection(1).data(point(2),point(1),:)));
% hold off;
% legend('real data','phantom');

%% Set up the phantom

mask_a0 = imfill(getLargestObj3D(s.a(1).data > 0.1),'holes');
a0 = s.a(1).data.*mask_a0;

% a0 = permute(a0,[2 3 1]);
% a0 = flip(a0,1);

%% define shale phantom parameters
N = 9; %number of projection sets to make with phi_op ranging from 90 degrees to 0 degrees.
numOfsegments = 32; %from 0 to pi
%Orientation distribution parameters
a0 = a0*0.7;
a2 = 1/5*0.7.*a0;
a4 = 1/48*0.7.*a0;
a6 = 1/96*0.7.*a0;

proj_data_size = [size(a0,1) size(a0,2)];
tomo_sz = size(a0);

savepath = [pwd '/'];
useID15Aangles = 1;
phi_op = deg2rad(80);%linspace(pi/2,0,N);
theta_op = deg2rad(90);
sigma = deg2rad(0);
%% Start loop here
for ii = 1:length(phi_op)
    
    phi_op_now = squeeze(squeeze(phi_op(ii)));
    
    %make distribution of angles
%     phi_op_3D = ones(size(a0)).*normrnd(phi_op_now,sigma,size(a0));
    phi_op_3D = ones(size(a0)).*phi_op_now; %.*normrnd(phi_op_now,sigma,size(a0));
    theta_op_3D = ones(size(a0)).*theta_op;
    phantom_3D = zeros(tomo_sz(1),tomo_sz(2),tomo_sz(3));
    
    phantom_3D_all_param = zeros([tomo_sz 5]);
    
    %set theta_op
    phantom_3D_all_param(:,:,:,1) = a0;
    phantom_3D_all_param(:,:,:,2) = a2;
    phantom_3D_all_param(:,:,:,3) = a4;
    phantom_3D_all_param(:,:,:,4) = a6;
    phantom_3D_all_param(:,:,:,5) = phi_op_3D;
    phantom_3D_all_param(:,:,:,6) = theta_op_3D;
    
    a0_lin = reshape(squeeze(phantom_3D_all_param(:,:,:,1)),[1 tomo_sz(1)*tomo_sz(2)*tomo_sz(3)]);
    a2_lin = reshape(squeeze(phantom_3D_all_param(:,:,:,2)),[1 tomo_sz(1)*tomo_sz(2)*tomo_sz(3)]);
    a4_lin = reshape(squeeze(phantom_3D_all_param(:,:,:,3)),[1 tomo_sz(1)*tomo_sz(2)*tomo_sz(3)]);
    a6_lin = reshape(squeeze(phantom_3D_all_param(:,:,:,4)),[1 tomo_sz(1)*tomo_sz(2)*tomo_sz(3)]);
    phi_op_lin = reshape(squeeze(phantom_3D_all_param(:,:,:,5)),[1 tomo_sz(1)*tomo_sz(2)*tomo_sz(3)]);
    theta_op_lin = reshape(squeeze(phantom_3D_all_param(:,:,:,6)),[1 tomo_sz(1)*tomo_sz(2)*tomo_sz(3)]);
    
    a = [a0_lin' a2_lin' a4_lin' a6_lin'];
    
    % Define projections angles
    if useID15Aangles
        [rot_x,rot_y] = get_ID15_SASTT_angles(); %using the projection angles in the Oct 2018 ID15A experiment.
    end
    
    % generate dummy projection to send into forward proj algorithm
    load('G:/PhD/Work/Archive/bone_xrdct/data/result_s3t2_reprocessed_09Dec2019_q17-40_all_again_add_id_rcoeff_0.105.mat');
%     load([pwd '/result_shale_q_38_q17-40_all_again_test_cap_reg_coeff_0.000.mat']);
    projection = struct();
    for jj = 1:length(rot_x)
        projection(jj).data = zeros(proj_data_size(1),proj_data_size(2),numOfsegments);
        projection(jj).rot_x = rot_x(jj);
        projection(jj).rot_y = rot_y(jj);
        projection(jj).Rot_exp = calcRotExp(rot_x(jj),rot_y(jj),0,0);
        projection(jj).window_mask = ones(proj_data_size(1),proj_data_size(2));
    end
    params.theta_0 = theta_op_lin;
    params.theta_init = theta_op_lin;
    params.phi_0 = phi_op_lin;
    params.phi_init = phi_op_lin;
    params.a_out1 = a;
    params.a = a;
    params.numOfvoxels = length(theta_op_lin);
    params.lines = tomo_sz(2);
    params.points_per_line = tomo_sz(1);
    params.regularization_angle = 0;
    params.numsegments = numOfsegments;
    params.phi_det = deg2rad(linspace(360/(4*numOfsegments),360-360/(4*numOfsegments),2*numOfsegments));
    params.phi_det = params.phi_det(1:32);
    clear proj_1_SH
    
    %Calculate the forward projection
    proj_SH = plot_projection_reg(params,projection(1:61));
    
    plotIt = 0;
    if plotIt
        point = [27 27];
        figure;
        plot(squeeze(proj_SH(20).projection(point(1),point(2),:)));
        axis equal;
    end
    
    clear projection_phantom
    for kk = 1:length(proj_SH)
        projection_phantom(kk).data = proj_SH(kk).projection;
        projection_phantom(kk).dx = 0;
        projection_phantom(kk).dy = 0;
        projection_phantom(kk).rot_x = proj_SH(kk).rotx;
        projection_phantom(kk).rot_y = proj_SH(kk).roty;
        projection_phantom(kk).Rot_exp = calcRotExp(proj_SH(kk).rotx,proj_SH(kk).roty,0,0);
        projection_phantom(kk).diode = squeeze(mean(projection_phantom(kk).data,3)).*0.001;
        projection_phantom(kk).par = getParForSASTT('s3t2'); %these are just dummy parameters that follow the data in the cSAXS pipeline
        projection_phantom(kk).integ.phi_det = linspace(360/(4*numOfsegments),360-360/(4*numOfsegments),2*numOfsegments);
        projection_phantom(kk).integ.phi_det = projection_phantom(kk).integ.phi_det(1:32);
        projection_phantom(kk).integ.norm_sum_avg = ones(numOfsegments,24);
        projection_phantom(kk).integ.norm_sum = ones(1233,numOfsegments);
    end

    % Save data
    savename = [savepath '/shale_phantom_test_20May2021.mat'];
%     savename = [savepath sprintf('/shale_phantom_28Aug2020_phi_op_mean_%0.2f_single_axis_aligned_ASTRA.mat',phi_op_now)];
    projection = projection_phantom;
    save(savename,'projection','phantom_3D_all_param');
    
end
