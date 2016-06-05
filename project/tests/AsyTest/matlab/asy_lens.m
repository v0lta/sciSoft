% Evaluate the asymptotic expansion for the generalised Jacobi polynomial 
% in the lens-shaped region with a certain normalisation. Derivative is only
% computed when assigned, default is monic normalisation without derivatives.
% Input
%   n,z          - Degree and point at which to evaluate
%   alpha,beta,h - Parts of the weight 
%   psi          - Anonymous function for analytic phase function psi(x)
%   nrT          - Number of terms. 1 gives leading order term
%   Dinf         - Limit of the Szego function
%   Uright,Uleft - Right and left U-matrices
%   [nor         - Normalisation: 'm' for monic(default), 'o' for orthonormal]
%   [dh,dpsi     - Derivatives of part of weight and phase function]
% Output
%   ipi          - Asymptotic expansion of the polynomial in the lens
%   [dipi        - Derivative of polynomial in the lens]
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created May 2014, last edit February 2015, shortened
%   verison.
function ipi = asy_lens(n,z,nrT,alpha,beta,Dinf,Uright,Uleft)

if (nrT == 1) || (numel(Uright) == 0), nmz = 1;	else
	nmz = sqrt(real(1+2*1i*Dinf^2*sum(reshape(Uright(2,1,1:(nrT-1),1) + ...
    Uleft(2,1,1:(nrT-1),1), [(nrT-1) 1])./(n+1).^(1:(nrT-1) )' ) ) );
end

RI = eye(2);
for k = 1:nrT-1
    for m = 1:round(k/2)
        RI = RI + (Uright(:,:,k,m)./(z-1).^m + Uleft(:,:,k,m)./(z+1).^m)/n.^k;
    end
end

ipi = [2.^(1/2-n) ./( (1+z).^(1/4+alpha/2).*(1-z).^(1/4+beta/2) ), 0]*RI*...
      [Dinf  *  cos((n+1/2)*acos(z) + (alpha+beta)/2 -  (alpha * pi)/2  -pi/4);...
      -1i/Dinf*cos((n-1/2)*acos(z) + (alpha+beta)/2 -  (alpha * pi)/2-pi/4)]*nmz;

end

