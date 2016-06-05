
%cd ./tests/AsyTest/matlab/

alpha = 0;
beta = 0;
maxP2 = 4;
maxOrder = 4;

Dinf = 2^(-alpha/2-beta/2);

[Uright,Uleft] = UExplicitToThree(alpha,beta,Dinf);

z = 0.2+0.5*1i;
n = 64;
nrT = 2;

ipi = zeros(maxP2,maxOrder);
for tn = 1:maxP2
    for index = 1:maxOrder
    ipi(tn,index) = asy_lens(tn^2,z,index,alpha,beta,Dinf,Uright,Uleft);
    end
end
disp(ipi)
exit
