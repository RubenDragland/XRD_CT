
%Test
x0 = dlarray([-1,2]);
tic
[fval,gradval] = dlfeval(@rosenbrock,x0);
toc
display(gradval);

x0 = [-1,2];
tic
[E, grad] = autograd_test(@rosenbrock, x0);
toc
display(grad);

% Parameters

mask_omega = ones(1,1,1);

a_coeffs = ones( 1,1,1,1,1) .* 69; 
Y_lms = ones(1,1,1,1,1) .*169 ; 

estimated_I = estimate_I(a_coeffs, Y_lms);
measured_I = ones(1,1,1,1) .* 19;
transparency = ones(1,1,1) .*0.5;

times = 10000;
% Time symbolic

tic
for ii = drange(1:times)
    sym_eps_grad = epsilon_coeff(mask_omega, mask_omega, estimated_I, measured_I, transparency, Y_lms, a_coeffs, Y_lms );
end
toc
display(sym_eps_grad);

% Time AD

tic
for ii = drange(1:times)
    [ad_eps, ad_grad] = matlab_autograd(@cost_function, a_coeffs, mask_omega, Y_lms, measured_I, transparency);
end
toc

display(ad_grad);


%Symbolic case

function [mask] = mask_term_gradient(mask_M, mask_omega)

mask = 4* mask_M .* mask_omega;
end

function [I] = intensity_term_gradient(estimated_I, measured_I, transparency)

I = ( estimated_I .^ 0.5 - (measured_I ./ transparency) .^ 0.5) ./ estimated_I .^0.5;
end

function [sum_SH] = sum_of_spherical_harmonics( a_coeffs, Y_lms)

sum_SH = sum(a_coeffs .* Y_lms, [4,5] ) ;
end

function [eps_q] = epsilon_coeff(mask_M, mask_omega, estimated_I, measured_I, transparency, Y_lm, a_coeffs, Y_lms )

eps_q = mask_term_gradient(mask_M, mask_omega) .* intensity_term_gradient(estimated_I, measured_I, transparency); % Elementwise multiplication
eps_q = sum( eps_q, [1,5]) ; % Matlab sum over axis? Remember 1-indexing in matlab...
eps_q = eps_q .* Y_lm(:,:); % Matlab indexing is unknown. 
eps_q = eps_q .* sum_of_spherical_harmonics(a_coeffs, Y_lms);
end

% Matlab autograd package
function [I_n] = estimate_I(a_coeffs, Y_lms)
I_n = sum(a_coeffs .* Y_lms, [4,5]); % Matlab sum over axis?
I_n = abs(I_n) .^2;
end

function [eps_q, d_eps_da] = cost_function(mask_omega, Y_lms, measured_I, transparency, a_coeffs)

estimated_I = estimate_I(a_coeffs, Y_lms);
eps_q = 2 * sum( mask_omega .* ( sqrt( estimated_I ) - sqrt( measured_I / transparency) ) .^2) ;
d_eps_da = dlgradient(eps_q, a_coeffs);
end

function[E, grad] =  matlab_autograd(func, diff_var, mask_omega, Y_lms, measured_I, transparency) %varargin)
x_diff = dlarray(diff_var);
%mask_omega = dlarray(mask_omega);
%Y_lms = dlarray(Y_lms);
%measured_I = dlarray(measured_I);
%transparency = dlarray(transparency);
[E,grad] = dlfeval(func, mask_omega, Y_lms, measured_I, transparency, x_diff); %varargin(:), x_diff); %Might need @ in front of function.
end

function[E, grad] = autograd_test(func, x)
x = dlarray(x);
[E, grad] = dlfeval(func, x);
end

% Test function
function [y,dydx] = rosenbrock(x)

y = 100*(x(2) - x(1).^2).^2 + (1 - x(1)).^2;
dydx = dlgradient(y,x);

end



