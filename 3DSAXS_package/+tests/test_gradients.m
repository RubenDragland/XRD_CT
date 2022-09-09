%% Test gradients
% This function should be improved by one that can create test data and
% check that the computed gradients are correct

clear
addpath utilities/
% Load data and optimization output

load('/das/work/units/csaxs/p17283/Bone_sample/analysis/SASTT/bone/projection_data/SASTT_bone_aligned_ASTRA.mat');
load('/das/work/units/csaxs/p17283/Bone_sample/analysis/SASTT/bone/SH/test_update-p-s-struct_20190529/optimization_output/result_bone_q10-20_angles_test_update-p-s-struct_20190529.mat');


if ~isfield(projection(1).integ,'theta_det')
    projection(1).integ.theta_det = pi/2;
end

%%%% Redefining optimization task %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test 1, gradient of first coefficient % 
p.find_orientation = false;
p.find_coefficients = true;
p.opt_coeff = [1 0 0];
opt = s.a(1).data;

indextoprobe = 269700;
value = s.a(1).data(indextoprobe);
delta = value*1e-3;

[E, grad, proj_out, Ereg] = optimization.SAXS_tomo_3D_err_metric(opt, p, s, projection );

s.a(1).data(indextoprobe) = value + delta;
opt = s.a(1).data;

[E2, grad2] = optimization.SAXS_tomo_3D_err_metric(opt, p, s, projection );
value2 = s.a(1).data(indextoprobe);
gradnum = (E2-E)/delta;
%
fprintf('*******************************************\n')
fprintf('  Test 1, gradient of first coefficient \n')
fprintf('Value of tensor 1 voxel %s = %f \n',indextoprobe,value)
fprintf('Value of tensor 2 voxel %s = %f \n',indextoprobe,value2)
fprintf('Value of analytic  gradient %s = %f \n',indextoprobe,grad(indextoprobe))
fprintf('Value of numerical gradient %s = %f \n',indextoprobe,gradnum)

figure(1)
plot(opt)
figure(2)
plot(grad)

%% Test 2, gradient of second coefficient % 
p.find_orientation = false;
p.find_coefficients = true;
p.opt_coeff = [1 1 0];
opt = s.a(1).data;
opt = [opt s.a(2).data];

indextoprobe = 929800;
value = opt(indextoprobe);
delta = value*1e-3;

[E, grad, proj_out, Ereg] = optimization.SAXS_tomo_3D_err_metric(opt, p, s, projection );

opt(indextoprobe) = value + delta;

[E2, grad2] = optimization.SAXS_tomo_3D_err_metric(opt, p, s, projection );
value2 = opt(indextoprobe);
gradnum = (E2-E)/delta;
%
fprintf('*******************************************\n')
fprintf('  Test 2, gradient of second coefficient \n')
fprintf('Value of tensor 1 voxel %s = %f \n',indextoprobe,value)
fprintf('Value of tensor 2 voxel %s = %f \n',indextoprobe,value2)
fprintf('Value of analytic  gradient %s = %f \n',indextoprobe,grad(indextoprobe))
fprintf('Value of numerical gradient %s = %f \n',indextoprobe,gradnum)

figure(1)
plot(opt)
figure(2)
plot(grad)
%% Test 3, gradient of third coefficient % 
p.find_orientation = false;
p.find_coefficients = true;
p.opt_coeff = [1 1 1];
opt = s.a(1).data;
opt = [opt s.a(2).data];
opt = [opt s.a(3).data];

indextoprobe = 1492000;
value = opt(indextoprobe);
delta = value*1e-3;

[E, grad, proj_out, Ereg] = optimization.SAXS_tomo_3D_err_metric(opt, p, s, projection );

opt(indextoprobe) = value + delta;

[E2, grad2] = optimization.SAXS_tomo_3D_err_metric(opt, p, s, projection );
value2 = opt(indextoprobe);
gradnum = (E2-E)/delta;
%
fprintf('*******************************************\n')
fprintf('  Test 3, gradient of third coefficient \n')
fprintf('Value of tensor 1 voxel %s = %f \n',indextoprobe,value)
fprintf('Value of tensor 2 voxel %s = %f \n',indextoprobe,value2)
fprintf('Value of analytic  gradient %s = %f \n',indextoprobe,grad(indextoprobe))
fprintf('Value of numerical gradient %s = %f \n',indextoprobe,gradnum)

figure(1)
plot(opt)
figure(2)
plot(grad)

%% Test 4, gradient of theta % 
p.find_orientation  = true;
p.find_coefficients = false;
p.opt_coeff = [0 0 0];
opt = s.theta.data;
opt = [opt s.phi.data];


indextoprobe = 352600+5404;
value = opt(indextoprobe);
delta = value*1e-3;

[E, grad, proj_out, Ereg] = optimization.SAXS_tomo_3D_err_metric(opt, p, s, projection );

opt(indextoprobe) = value + delta;

[E2, grad2] = optimization.SAXS_tomo_3D_err_metric(opt, p, s, projection );
value2 = opt(indextoprobe);
gradnum = (E2-E)/delta;
%
fprintf('*******************************************\n')
fprintf('  Test 4, gradient of theta \n')
fprintf('Value of tensor 1 voxel %s = %f \n',indextoprobe,value)
fprintf('Value of tensor 2 voxel %s = %f \n',indextoprobe,value2)
fprintf('Value of analytic  gradient %s = %f \n',indextoprobe,grad(indextoprobe))
fprintf('Value of numerical gradient %s = %f \n',indextoprobe,gradnum)

figure(1)
plot(opt)
figure(2)
plot(grad)

%% Test 4, gradient of phi % 
p.find_orientation  = true;
p.find_coefficients = false;
p.opt_coeff = [0 0 0];
opt = s.theta.data;
opt = [opt s.phi.data];


indextoprobe = 902700+36;
value = opt(indextoprobe);
delta = value*1e-3;

[E, grad, proj_out, Ereg] = optimization.SAXS_tomo_3D_err_metric(opt, p, s, projection );

opt(indextoprobe) = value + delta;

[E2, grad2] = optimization.SAXS_tomo_3D_err_metric(opt, p, s, projection );
value2 = opt(indextoprobe);
gradnum = (E2-E)/delta;
%
fprintf('*******************************************\n')
fprintf('  Test 5, gradient of phi \n')
fprintf('Value of tensor 1 voxel %s = %f \n',indextoprobe,value)
fprintf('Value of tensor 2 voxel %s = %f \n',indextoprobe,value2)
fprintf('Value of analytic  gradient %s = %f \n',indextoprobe,grad(indextoprobe))
fprintf('Value of numerical gradient %s = %f \n',indextoprobe,gradnum)

figure(1)
plot(opt)
figure(2)
plot(grad)
