clear all;

%define a coarse grid.
x = load('x.in');
%compute the sin on it.
f = load('fun.in');

[x,I] = sort(x);
f = f(I);

%define a finer grid.
des = load('des.in');
des = sort(des);


%TODO: rescale to input to -1,1
%compute the transformation parameters:
b = 2/(max(x) - min(x));
a = (-1/b) -min(x);

x = (x + a).*b;

alpha = 0;
beta = 0;
M = 4;
maxP2 = M;
maxOrder = M;

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
orthPol = real(ipi');

%find the p matrix:
%figure(1)
for i = 1:M
    P(:,i) = polyval(orthPol(i,:),x);
    %plot(x,P(:,i)); hold on;
end
%hold off;

c = pinv(P) * f';

des = (des + a).*b;

for i = 1:M
    P2(:,i) = polyval(orthPol(i,:),des);
end
orthRep = P2*c;

%back transformation
x = x/b - a;
des = des/b -a;


sol = tanh(des) + besseli(3,des);

figure(2)
plot(des,orthRep,'*');
hold on;
plot(des,sol,'*');
hold off;
%plot(c)
