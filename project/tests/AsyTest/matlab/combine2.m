
%cd ./tests/AsyTest/matlab/

alpha = 0.2;
beta = -0.2;
pointNo = 4;
maxOrder = 4;

Dinf = 2^(-alpha/2-beta/2);

[Uright,Uleft] = UExplicitToThree(alpha,beta,Dinf);


z = linspace(-0.9,0.9,pointNo);
%z = sqrt(eps) - 1 + (z - min(z)) / (max(z) - min(z)) * (2 - 2 *sqrt(eps));
n = 6;

orthPol = zeros(pointNo,maxOrder);
for i = 1:pointNo
    for j = 1:maxOrder
        orthPol(j,i) = asy_lens(n,z(i),j,alpha,beta,Dinf,Uright,Uleft);
    end
end
disp(orthPol)
exit
