% Get the U-matrices to construct the asymptotic expansion of R using the
% procedure with the convolutions as explained in the paper. Optionally,
% specify method or number of matrices to compute or get the Q-s.
% Input
%   alpha, beta  - Parts of the weight w(x) = (1-x).^alpha.*(1+x).^beta.*h(x)
%   Dinf         - Limit of the Szego function
%   c,d          - Coefficients needed for the U-matrices
%   maxOrder     - The maximum order of the error
% Output
%   Uright       - Coefficient matrices of R_k(z) for (z-1)^(-m)
%   Uleft        - Coefficient matrices of R_k(z) for (z+1)^(-m)
%   [Qright      - Coefficient matrices of R_k^{right}(z) for (z-1)^n]
%   [Qleft       - Coefficient matrices of R_k^{left}(z) for (z+1)^n]
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created January 2014, last edit February 2015,
%   simplified version
function [Uright,Uleft,Qright,Qleft] = UQ(alpha,beta,Dinf,c,d,maxOrder)

Uright = zeros(2,2,maxOrder-1,ceil( (maxOrder-1)/2) );
Uleft = zeros(2,2,maxOrder-1,ceil( (maxOrder-1)/2));
% About half of these tensors will not be used

Qright = zeros(2,2,maxOrder-1,ceil( (maxOrder-1)/2) );
Qleft = zeros(2,2,maxOrder-1,ceil( (maxOrder-1)/2));


poch = @(x,n) prod(x+(0:(n-1) ) ); % 'pochhammer' needs symbolic toolbox, 
% could also use separate 1-line script or global function

Wr = WV(alpha,beta,Dinf,c,maxOrder,1,1);
Wl = WV(beta,alpha,Dinf,d,maxOrder,-1,1);

for k = 1:(maxOrder-1)
		% Actually, for m = (-round(k/2)):1
		% but Matlab does not allow negative indices so shift by round(k/2)+1
		for m = 1:round(k/2)
			Uright(:,:,k,m) = Wr(:,:,k,round(k/2)+1-m);
			Uleft(:,:,k,m) = Wl(:,:,k,round(k/2)+1-m);
			for j=1:(k-1)
				for l = max(m-ceil(j/2),1):ceil((k-j)/2)
					Uright(:,:,k,m) = Uright(:,:,k,m) + ...
						Uright(:,:,k-j,l)*Wr(:,:,j,round(j/2)+1+l-m);
					Uleft(:,:,k,m) = Uleft(:,:,k,m) + ...
						Uleft(:,:,k-j,l)*Wl(:,:,j,round(j/2)+1+l-m);
				end
            end
        for j=1:(k-1)
			for n = 0:(ceil(j/2)-m)
				for i = 1:ceil((k-j)/2)
						Uright(:,:,k,m) = Uright(:,:,k,m) + ...
							poch(1-i-n,n)/2^(i+n)/factorial(n)*...
							Uleft(:,:,k-j,i)*Wr(:,:,j,round(j/2)+1-n-m);
						Uleft(:,:,k,m) = Uleft(:,:,k,m) + ...
							poch(1-i-n,n)/(-2)^(i+n)/factorial(n)*...
							Uright(:,:,k-j,i)*Wl(:,:,j,round(j/2)+1-n-m);
                end
			end
		end
	end
end

end
