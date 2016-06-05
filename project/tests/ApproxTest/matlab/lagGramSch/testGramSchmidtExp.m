clear all;

%define a coarse grid.
x = -pi:0.1:pi;
%compute the sin on it.
f = exp(x);


%define a finer grid.
des = [-pi/16, 0.0 , pi/16];



%TODO: rescale to input to -1,1
%compute the transformation parameters:
b = 2/(max(x) - min(x));
a = (-1/b) -min(x);

x = (x + a).*b;

alpha = 0;
beta = -0;

M = 3;
monPoly = flipud(eye(M));
orthPol = gramSchmidt( alpha,beta,monPoly );

%find the p matrix:
figure(1)
for i = 1:M
    P(:,i) = polyval(orthPol(i,:),x);
    plot(x,P(:,i)); hold on;
end
hold off;

c = pinv(P) * f';

des = (des + a).*b;

for i = 1:M
    P2(:,i) = polyval(orthPol(i,:),des);
end
orthRep = P2*c;

%back transformation
x = x/b - a;
des = des/b -a;

sol = load('solF.test');
sol = sort(sol);

figure(2)
plot(des,orthRep,'*');
hold on;
plot(des,exp(des),'*');
plot(x,exp(x));
hold off;
%plot(c)
