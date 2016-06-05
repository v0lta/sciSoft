function [ orthpol ] = gramSchmidt( alpha,beta,orthpol )
%compute a monomial orthonormal polynomial basis.

omega = @(x) (1-x).^alpha .* (1+x).^beta;
sze = size(orthpol);
rows = sze(1);


for i=1:rows
    
    %normalize:
    normFun = @(x) omega(x) .* polyval(orthpol(i,:),x).^2;
    orthpol(i,:) = orthpol(i,:) ./ sqrt(integral(normFun,-1,1));
        
    %orthogonalize
    for j=1:i
        if (i ~= j)
            numFun = @(x) omega(x).*polyval(orthpol(i,:),x).*polyval(orthpol(j,:),x);
            denomFun = @(x) omega(x).*polyval(orthpol(i,:),x).^2;
            
            parPart = real(sqrt(integral(numFun,-1,1)/integral(denomFun,-1,1)));
            orthpol(i,:) = orthpol(i,:) - parPart .* orthpol(j,:);
        end
    end 
end 


end

