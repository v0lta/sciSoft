clear all;
format long

%define a coarse grid.
x = -1*pi:0.1:1*pi;
f = @(x) sin(x);
%f = @(x) cos(x);

des = -1*pi:0.1:1*pi;
alpha = 0;
beta = 0;


N = 50;
ab = r_jacobi(N,alpha,beta);
xw = gauss(N,ab);


M = 5;
monPoly = flipud(eye(M));
orthPol = gramSchmidt( alpha,beta,monPoly );

omega = @(x) (1-x).^alpha .* (1+x).^beta;

for j = 1:M
    intFun = @(x) omega(x) .* polyval(orthPol(j,:),x) .* f(x);
    sum = 0;
    for i= 1:N
        sum = sum + intFun(xw(i,1))*xw(i,2);
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
plot(x,f(x))
