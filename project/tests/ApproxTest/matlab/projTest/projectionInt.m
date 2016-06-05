clear all;
format long

epsilon = 2.2204460492503131E-016;
%define a coarse grid.
x = -2:0.1:2;
f = @(x) sin(x);
%f = @(x) cos(x);

des = -2:0.01:2;
alpha = 0;
beta = 0;

M = 5;
orthPol = favard(M,alpha,beta);

omega = @(x) (1-x).^alpha .* (1+x).^beta;

epsilon = 2.2204460492503131E-016;
x = sqrt(epsilon) - 1 + (x - min(x))/ (max(x) - min(x)) * (2 - 2*sqrt(epsilon)) 

for j = 1:M
    intFun = @(x) omega(x) .* polyval(orthPol(j,:),x) .* f(x);
    sumVar = integral(intFun,-1,1);
    c(j) = sumVar;
end 

%find the p matrix:
figure(1)
for i = 1:M
    P(:,i) = polyval(orthPol(i,:),des);
    %plot(des,P(:,i)); hold on;
end
%hold off;

orthRep = P * c';

plot(des,orthRep,'*')
hold on;
plot(x,f(x))
