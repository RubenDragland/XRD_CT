
function [E, grad, proj_out, Ereg] = SAXS_EXPSIN_3D_err_AD(opt_inputs, p, s, projection)

if isempty(opt_inputs)
    skip_optimization = true;
else
    skip_optimization = false;
end

% RSD: Do not bother with parallel atm. 
% if isempty(gcp('nocreate'))
%     parpool('local', feature('numcores'));
% end
find_grad = (nargout >= 2);

if skip_optimization
    warning('Skipping optimization')
    grad = [];
    find_grad = false;
end


phi_det = p.phi_det;                            % Measurement azimuthal angles (detector)
theta_det = ones(size(phi_det)).*projection(1).integ.theta_det; 

nx = p.nx;                            
ny = p.ny;
nz = p.nz;
numOfsegments = p.numsegments;          % number of segements the saxs data is averaged (default is 16)
numOfvoxels = p.numOfvoxels;            % dimension of tomogram
mask3D = s.mask3D;
numOfpixels = numel(projection)*nx*ny;

%RSD: Prepare theta
if ~skip_optimization
    theta_struct = opt_inputs(1:numOfvoxels);
    phi_struct = opt_inputs(numOfvoxels+1:2*numOfvoxels);
    theta_struct = reshape(theta_struct,[p.ny,p.nx,p.nz]);
    phi_struct   = reshape(phi_struct  ,[p.ny,p.nx,p.nz]); 
    opt_inputs = opt_inputs(2*numOfvoxels+1:end);

    if p.avoid_wrapping
        phi_struct = mod(phi_struct, 2*pi);
        theta_struct = mod(theta_struct, 2*pi);
    end
    
    A = [];
    B = [];
    A = [A; opt_inputs(1:numOfvoxels)];
    opt_inputs = opt_inputs(numOfvoxels+1:end);
    B = [B; opt_inputs(1:numOfvoxels)];
else
    theta_struct = s.theta.data;
    phi_struct = s.phi.data;
    
    A = s.A.data;
    B = s.B.data;
end

A = reshape(A, numOfvoxels, 1);
B = reshape(B, numOfvoxels, 1);

% grad_a = zeros(ny, nx, nz); %RSD: Bit unsure on shape. 
% 
% grad_theta_struct = zeros(ny, nx, nz);
% grad_phi_struct   = zeros(ny, nx, nz);

E = 0;

% q (unit vector), 3D coordinates of the measurements in reciprocal space and in beamline coordinates. 
% Size 3x8 [xyz x numOfsegments]. The third dimension is actually only needed to account for Ewald sphere curvature
% correction which is included in theta_det
unit_q_beamline = [sin(theta_det).*cos(phi_det); sin(theta_det).*sin(phi_det);cos(theta_det)]; 

N = [ny nx nz]; % Size of tomogram  
x = (1:N(2)) - ceil(N(2)/2); %RSD: Is everything fine with this?
y = (1:N(1)) - ceil(N(1)/2);
z = (1:N(3)) - ceil(N(3)/2);
[X, Y, Z] = meshgrid(x, y, z); % (x,y,z) (2,1,3)

P = py.sys.path; %RSD: Fix this path to not be hard-coded. 
search_path = [p.parent '/Autodiff_package'];
if count(P,search_path) == 0
    insert(P,int32(0),search_path);
end

py.importlib.import_module("batch_AD");

theta_struct_it = py.numpy.array(theta_struct);
phi_struct_it = py.numpy.array(phi_struct);
A_it = py.numpy.array(A);
B_it = py.numpy.array(B);
Xi = py.numpy.array(X);
Yi = py.numpy.array(Y);
Zi = py.numpy.array(Z);
unit_q_beamline_i = py.numpy.array(unit_q_beamline);

batch_results = py.batch_AD.EXPSIN_main( theta_struct_it, phi_struct_it, A_it, B_it, unit_q_beamline_i, p, Xi, Yi, Zi, ny, nx, nz, numOfsegments, numOfpixels, numOfvoxels, find_grad);

E = batch_results{1};

if ~isempty( double(batch_results{2}) )
    AD_grad_A = double(batch_results{2});
    AD_grad_B = double(batch_results{3});
    grad_A = reshape( permute(AD_grad_A, [3,2,1]) , ny, nx, nz) ; 
    grad_B = reshape( permute(AD_grad_B, [3,2,1]), ny,nx,nz);
end
if ~isempty( double(batch_results{4}) )
    AD_grad_theta = double( batch_results{4});
    AD_grad_phi = double (batch_results{5} );
    grad_theta_struct = reshape( permute(AD_grad_theta, [3,2,1]), ny, nx, nz); 
    grad_phi_struct = reshape( permute(AD_grad_phi, [3,2,1]), ny,nx,nz);
end

if find_grad
% gradient output vector
    grad = [];
    global Err_hist
    Err_hist = [Err_hist, E];

    grad_theta_struct = reshape(grad_theta_struct, 1, []);
    grad_phi_struct = reshape(grad_phi_struct, 1, []);
    grad_A = reshape(grad_A, 1, []);
    grad_B = reshape(grad_B, 1, []);
    grad = [grad, grad_theta_struct, grad_phi_struct, grad_A, grad_B];
end

if find_grad
    iteration = length(Err_hist);
    fprintf('*************************************************** \nIteration %d \n', iteration);
else
    fprintf('Line Search \n')
end
fprintf('          data_error = %f \n', E);


end