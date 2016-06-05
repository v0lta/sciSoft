clear all;
format long

epsilon = 2.2204460492503131E-016;
%define a coarse grid.
x = -pi:0.1:pi;
f = exp(x);


des = x;
alpha = 0;
beta = 0;


N = length(f);
ab = r_jacobi(N,alpha,beta);
xw = gauss(N,ab);

epsilon = 2.2204460492503131E-016;
x = sqrt(epsilon) - 1 + (x - min(x))/ (max(x) - min(x)) * (2 - 2*sqrt(epsilon));
des = sqrt(epsilon) - 1 + (des - min(des))/ (max(des) - min(des)) * (2 - 2*sqrt(epsilon));


M = 9;
monPoly = flipud(eye(M));
orthPol = gramSchmidt( alpha,beta,monPoly );

omega = @(x) (1-x).^alpha .* (1+x).^beta;

for j = 1:M
    sum = 0;
    for i= 1:N
        intFun = omega(x(i)) .* polyval(orthPol(j,:),x(i)) .* f(i);
        sum = sum + intFun*xw(i,2);
    end
    c(j) = sum;
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
plot(x,f)
