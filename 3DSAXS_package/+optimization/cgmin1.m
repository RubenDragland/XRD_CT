function [x] = cgmin1(func,x,itmax,ftol,xtol,varargin)
% conjugate-gradient optimization routine
% NOTE: linesearch subroutines do not use the gradient
%
% [x] = cgmin1(func,x,itmax,ftol,xtol,varargin)
%
%     func = string name of objective function which returns both the
%            objective function value and the gradient
%        x = input as initial starting point and output as final point
%    itmax = maximum number of iterations (empty for default = 50)
%     ftol = relative function tolerance (empty for default = 1e-3)
%     xtol = absolute solution tolerance (empty for default = 1e-3)
% varargin = extra variables required by objective function
%
% DISCLAIMER: This code is not intended for distribution. I have many
% versions of this code and am constantly revising it. I believe this
% version is working properly. However, I will not vouch for the code.
% Anyone using the code for thesis research has a responsibility to go
% through the code line-by-line and read relevant references to understand
% the code completely. In my opinion, you have two options if you want to
% publish results obtained with the code: (i) go through the code line-by-
% line and read relevent references to understand how the code works and make
% sure it is working properly for your application, or (ii) I can sit down
% with you an go through this code and the additional code that you have
% written to go along with it and make sure it is working properly. Option
% (i) is preferred, and I ask that you do NOT acknowledge me in print (first,
% it would be more appropriate for you to reference "Numerical Recipes",
% and second, I prefer not to be named in a paper with which I do not have
% detailed knowledge). If you decide to go with option (ii), I would expect
% to learn the details of your research and be included in the author list.
%
% Sam Thurman, May 9, 2005
if isempty(itmax), itmax = 50; end
if isempty(ftol), ftol = 1e-3; end
if isempty(xtol), xtol = 1e-3; end
% loop
flg = 0; % use steepest descent for first iteration
step = 0; % to guess at initial steplength
for it = 1:itmax
   % function evaluation
   [f,grad] = feval(func,x,varargin{:});
%    disp(f)
   % check for feasibility
   if isinf(f), error('encountered an infeasible solution'), end
   if norm(grad(:))==0, return, end % done if gradient is zero (unlikely)
   % pick search direction
   if (flg==1) & (rem(it,25)~=0) % linesearch found a minimum -> use cg equations
       gg = g(:)'*g(:);
%        dgg = grad(:)'*grad(:); % this statement for Fletcher-Reeves
       dgg = (grad(:)+g(:))'*grad(:);  % this statement for Polak-Ribiere
       ga = dgg/gg;
       g = -grad;
       h = g+ga*h;
       dx = h/norm(h(:));
       df = grad(:)'*dx(:);
   end
   if (flg==0) | (rem(it,25)==0) | (df>0) % revert to steepest decent
       g = -grad;
       h = g;
       dx = h/norm(h(:));
       df = grad(:)'*dx(:);
   end
   % initial steplength guess
   if step == 0
       step = max(0.001,min([1,2*abs(f/(grad(:)'*dx(:)))])); % same as fminusub.m (line 124) in optim toolbox
   else % oterwise use previous steplength
        step = step/10;
   end
   % linesearch
   [x,fvalue,step,flg] = linesearch(func,x,f,df,dx,step,varargin{:});
   % test for convergence
   if (2*abs(f-fvalue)<=ftol*(abs(f)+abs(fvalue)+ftol)) & (step*norm(dx(:))<=xtol) & (it~=1) % normal return
       return
   end
end
disp('maximum number of iterations exceeded')
return

function [x,f,a,flg] = linesearch(func,x0,f0,df0,dx,a,varargin)
% linesearch routine (does not use the gradient)
%
% [x,f,a,flg] = linesearch(func,x0,f0,df0,dx,a,varargin)
%
% func = string name of objective function
% x0 = starting point of search
% f0 = objective function at x0
% df0 = derivative of objective function along dx at x0
% dx = direction of linesearch
% a = steplength input as guess output as taken
% x = final point of search
% f = objective function at final point
% flg = indicates how step was determined (0 for Armijo step, 1 for
%        bracketing and refining a minimum)

% check if descent direction
if df0>=0
   warning('linesearch called w/o descent direction')
   x = x0; f = f0; fcount = [0,0]; flg = 0;
   return
end
% first point
a1 = 0; f1 = f0;
% try initial steplength
f = feval(func,x0+a*dx,varargin{:}); % no gradient returned
% keep track of old values & hopefully bracket a minimum
a2 = a; f2 = f;
% make sure initial step is feasible
while isinf(f)
   a2 = a;
   a = 0.25*a;
   f = feval(func,x0+a*dx,varargin{:});
end
% decide what to do next based on Armijo condition
b = 0; % parameter in Armijo condition (>=0)
% if f does not satisfy Armijo condition (steplength may be too
% large) -> decrease steplength until Armijo is satisfied
if f>f0+a*b*df0
   while f>f0+a*b*df0
       if isinf(f2)|(f<=f2)
           a2 = a; f2 = f;
       end
       a = 0.25*a; % decrease steplength
       f = feval(func,x0+a*dx,varargin{:});
   end
end
% arrive here with a1=a0=0 and f<f1 (at least)
tmp = 1;
while f2<=f % try doubling a2
   a2 = 2*a2;
   f2 = feval(func,x0+a2*dx,varargin{:});
   if f2<f
       a1 = a; f1 = f;
       a = a2; f = f2;
   end
   if tmp==5; break; else tmp = tmp+1; end
end
% case where we should have a bracket, but f2 is infinite
while isinf(f2)
%    disp('should have a bracket but f2 is infinite')
   u = a+0.25*(a2-a); % point between a and a2
   fu = feval(func,x0+u*dx,varargin{:});
   if fu<f
       a1 = a; f1 = f;
       a = u;  f = fu;
   else
       a2 = u; f2 = fu;
   end
end
% last steps
if (f<f1)&(f<f2)&isfinite(f2) % bracketing successful -> refine minimum
   [a,f] = brent(func,x0,dx,a1,a,a2,f1,f,f2,varargin{:});
   flg = 1; % use conjugate gradient next loop
else % bracketing unsuccessful -> stop
   flg = 0; % use steepest descent next loop
end
x = x0+a*dx;
return

function [a,f] = brent(func,x0,dx,a1,a,a2,f1,f,f2,varargin)
% one-dimensional minimization by parabolic interpolation & golden
% section (does not use the gradient)
%
% [a,f] = brent(func,x0,dx,a1,a,a2,f1,f,f2,varargin)
%
% func = string name of objective function
% x0 = starting point of linesearch
% dx = direction of linesearch
% a1,a,a2 = bracketing triplet of steplengths (a1<a<a2)
% f1,f,f2 = objective function at steplengths a1, a, & a2
% a = output as final steplength
% f = output as objective function at final steplength

gold = 0.3819660; % golden ratio
itmax = 5;
tol = 0.5;
% check order of bracket
if (a1>a)|(a2<a), error('brent called with bracket in wrong order'), end
% initialize
v = a;fv = f; % middle point on step before last
w = a;fw = f; % middle point on last step
e = 0; % distance moved on step before last
% iterations
for it = 1:itmax
   am = 0.5*(a1+a2);
   tol1 = tol*abs(a)+eps;
   tol2 = 2*tol1;
   % test for convergence
   if abs(a-am)<=(tol2-0.5*(a2-a1)), return, end
   % choose next point
   if abs(e)>tol1 % construct a trial parabolic fit
       r = (a-w)*(f-fv);
       q = (a-v)*(f-fw);
       p = (a-v)*q-(a-w)*r;
       q = 2*(q-r);
       if q>0, p = -p; end
       q = abs(q);
       etemp = e;
       e = d;
       % check acceptability of parabolic fit
       ok = ~(abs(p)>=abs(0.5*q*etemp) | p<=q*(a1-a) | p>=q*(a2-a));
       if ok % take parabolic step
           d = p/q;
           u = a+d;
           if (u-a1)<tol2 | (a2-u)<tol2, d = sign(am-a)*tol1; end
       else % take golden section step
           if a>=am
               e = a1-a;
           else
               e = a2-a;
           end
           d = gold*e;
       end
   else % take golden section step
       if a>=am
           e = a1-a;
       else
           e = a2-a;
       end
       d = gold*e;
   end
   % arrive here with d computed either from
   % parabolic fit or else from golden section
   if abs(d)>=tol1
       u = a+d;
   else
       u = a+sign(d)*tol1;
   end
   fu = feval(func,x0+u*dx,varargin{:}); % one function evaluation per iteration
   if fu<=f
       if u>=a
           a1 = a;
       else
           a2 = a;
       end
       v = w; fv = fw;
       w = a; fw = f;
       a = u; f = fu;
   else
       if u<a
           a1 = u;
       else
           a2 = u;
       end
       if fu<=fw | w==a
           v = w; fv = fw;
           w = u; fw = fu;
       elseif fu<=fv | v==a | v==w
           v = u; fv = fu;
       end
   end
end
disp('exceeded maximum number of iterations')
return

