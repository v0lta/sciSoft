% Compute the contour integrals c_n, d_n, D_\infty, psi(z) [, its derivative and its contour integral].
% Input
%   alpha, beta - Parts of the weight w(x) = (1-x).^alpha.*(1+x).^beta.*h(x)
%   h           - Anonymous function for analytic function h(x)
%   nrT         - Number of terms, where 1 gives leading order term
%   [rho        - Berstein ellipse parameter: small to not encircle where log(h) is 
% 				  analytic and large to encircle the point x where to evaluate psi(x)]
%   [M          - Number of points in trapezoidal rules]
% Output
%   c, d        - The c- and d-coefficients, size(c&d) == ceil(nrT/2)
%   Dinf        - Limit of the Szego function
%   psi         - Anonymous function for the phase function
%   [dpsi       - Derivative of the phase function]
%   [contpsi	- Contour integral in psi, contpsi(x,p) has (zeta-x)^p]
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created October 2013, last edit February 2015
function [c, d, Dinf, psi, dpsi,contpsi] = contour_integrals(alpha,beta,h,nrT,rho,M)
c = zeros(ceil(nrT/2),1);
d = zeros(ceil(nrT/2),1);

ra = 2*rand-1; 
mm = 2; cc = 7; 
% Checking for known exact solutions by comparing value at one random z =
% ra, could also do if sum(chebfun(@(x) abs(h(x) -ones(size(x) ) ), domain(xc) ) ) == 0
if h(ra) == 1
    display('Using exact result for non-generalised weight in contour_integrals')
    Dinf = 2^(-alpha/2-beta/2);
    psi = @(x) 1/2.*( (alpha + beta).*acos(x)-alpha*pi);
    dpsi = @(x) 1/2.*(alpha+beta).*(-1./(sqrt(1-x).*sqrt(1+x) ) );
	contpsi = @(x,p) 0;
elseif abs(h(ra) - exp(-cc*ra.^(2*mm)) ) <100*eps
    display('Using exact result for Gibbs weight in contour_integrals');
	binom = @(x,n) prod(x:-1:(x-n+1) )/prod(1:n);
	% Could go until length(c)-1, but need exact result for psi, so:
    for n = 0:(2*mm-1)
		c(n+1) = 0;
        for j = 0:floor( (2*mm-n-1)/2)
            c(n+1) = c(n+1) -cc*binom(j-1/2,j)*binom(2*mm-1-2*j,2*mm-n-1-2*j);
        end
        d(n+1) = (-1)^(n+1)*c(n+1);
    end
    Dinf = 2^(-alpha/2-beta/2)*exp(-cc/2*binom(mm-1/2,mm));
    degs = (0:(length(c)-1) ).';
    psi = @(x) 1/2.*( (alpha + beta).*acos(x)-alpha*pi) +sqrt(1-x).*sqrt(1+x)/2.*sum( c.*(x-1).^(degs) );
    dpsi = @(x) 1/2.*(alpha+beta).*(-1./(sqrt(1-x).*sqrt(1+x) ) ) + ...
        (-x)./sqrt(1-x)./sqrt(1+x)/2.*sum(c.*(x-1).^(degs) ) + ...
        sqrt(1-x).*sqrt(1+x)/2.*sum(degs(2:end).*c(2:end).*(x-1).^(degs(2:end)-1) );
	contpsi = @(x,p) 2i*pi*sum( (degs(p:end)).^(p-1).*c(p:end).*(x-1).^(degs(p:end)-p+1) );
else % Trapezoidal rules
	if abs(h(ra) - 1./sqrt(ra+3)) < 100*eps
		display('Using trapezoidal rules for Fourier weight in contour_integrals');
		rho = 4; M = 80;
	elseif ~exist('rho','var') % Take rho smaller to reduce risk of encircling pole log h:
		rho = 1.25; % > 1 so positive direction
	end
	if ~exist('M','var') 
		M = 400; % for when we pass closely by pole
	end
    for nr = 1:(ceil(nrT/2)-1) % nr = n+1 in c_n = c(nr)
        c(nr) = trap_rule(@(z) log(h(z ) )./( sqrt(z -1).*sqrt(z +1).*(z -1).^nr ), rho,M);
        d(nr) = trap_rule(@(z) log(h(z ) )./( sqrt(z -1).*sqrt(z +1).*(z +1).^nr ), rho,M);
    end
    Dinf = 2.^( -(alpha+beta)./2).*exp(real(trap_rule(@(z) log(h(z) )./sqrt(z -1)./sqrt(z +1), rho,M)/2) );
    
    psi = @(x) 1/2.*( (alpha + beta).*acos(x)-alpha*pi) + sqrt(1-x).*sqrt(1+x) ...
        .*trap_rule(@(z) log(h(z) )./sqrt(z-1)./sqrt(z+1)./(z-x), rho, M)/2;
    % trap_rule already devides by 2i*pi
    dpsi = @(x) 1/2.*(alpha + beta).*(-1)./(sqrt(1-x).*sqrt(1 +x)) + ...
            (-x)/(sqrt(1-x).*sqrt(1+x)*2)*trap_rule( @(z) log(h(z) )./sqrt(z-1)./sqrt(z+1)./(z-x),rho,M) +...
            sqrt(1-x).*sqrt(1+x)/2*trap_rule(@(z) log(h(z) )./sqrt(z-1)./sqrt(z+1)./(z-x).^2,rho,M);
	contpsi = @(x,p) 2i*pi*trap_rule(@(z) log(h(z) )./sqrt(z-1)./sqrt(z+1)./(z-x).^p, rho, M);
end
