
%cd ./tests/AsyTest/matlab/

alpha = 0;
beta = 0;
pointNo = 4;
maxOrder = 4;

Dinf = 2^(-alpha/2-beta/2);

[Uright,Uleft] = UExplicitToThree(alpha,beta,Dinf);


z = linspace(-0.5,0.5,pointNo);
%z = sqrt(eps) - 1 + (z - min(z)) / (max(z) - min(z)) * (2 - 2 *sqrt(eps));
n = 200;

orthPol = zeros(pointNo,maxOrder);
for i = 1:pointNo
    for j = 1:maxOrder
        orthPol(i,j) = asy_lens(i^2,z(i),j,alpha,beta,Dinf,Uright,Uleft);
    end
end
disp(orthPol)
