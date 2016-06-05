% Compute the contour integral of the given function using the trapezoidal rule around the interval.
% Input
%   f       - the integrand
%   rho     - parameter for the Bernstein ellipse
%   M       - number of terms
% Output
%   cont    - the contour integral
% About
%   Author  - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History - Created October 2013, last edit February 2015
function cont = trap_rule(f,rho,M)
ths = linspace(0,2*pi,M+1); % = theta_k, k = 0 to M-1
ths(end) = [];
hh = ths(2)-ths(1);
zs = rho/2*exp(1i*ths) + 1/rho/2*exp(-1i*ths);
cont = hh*sum(f(zs).*(rho*1i/2*exp(1i*ths) -1i/rho/2*exp(-1i*ths) ) )/(2*pi*1i );
if rho < 1 % Bernstein in negative direction but need positive direction
    cont = -cont;
end
% Below an ellipse, but could also need other parametrisations, for example
% when evaluating psi(x) for x further from the interval than the point
% where log(h) is not analytic
if 0
    r1 = 1.8;
    r2 = 0.8;
    zs = r1*cos(ths) +r2*1i*sin(ths);
    cont = hh*sum(f(zs).*(-r1*sin(ths) +r2*1i*cos(ths) ) )/(2*pi*1i );
end