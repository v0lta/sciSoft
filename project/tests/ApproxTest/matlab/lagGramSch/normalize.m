function [ orthpol ] = normalize( alpha,beta,orthpol )
%compute a monomial orthonormal polynomial basis.

omega = @(x) (1-x).^alpha .* (1+x).^beta;
sze = size(orthpol);
rows = sze(1);

for i=2:rows
    %normalize:
    normFun = @(x) omega(x) .* polyval(orthpol(i-1,:),x).^2;
    orthpol(i-1,:) = orthpol(i-1,:) ./ sqrt(integral(normFun,-1,1));
    
    %orthoginalize:
    %for j=1:(i-1)
    %    numFun = @(x) omega(x).*polyval(orthpol(i,:),x).*polyval(orthpol(j,:),x);
    %    denomFun = @(x) omega(x).*polyval(orthpol(j,:),x).^2;
    %    
    %    parPart = integral(numFun,-1,1)/integral(denomFun,-1,1) * orthpol(j,:);
    %    orthpol(i,:) = orthpol(i,:) - parPart;
    %end
end 

%normalize:
normFun = @(x) omega(x) .* polyval(orthpol(end,:),x).^2;
orthpol(end,:) = orthpol(end,:) ./ sqrt(integral(normFun,-1,1));
end
