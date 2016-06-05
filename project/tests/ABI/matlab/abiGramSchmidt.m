clear all;

%define a coarse grid.
f = load('abiclose.in');
f = fliplr(f);
x = linspace(-1,1,length(f));

des = [linspace(-1,1,length(f)*2) 1.1 1.2 1.3];


%TODO: rescale to input to -1,1
%compute the transformation parameters:
b = 2/(max(x) - min(x));
a = (-1/b) -min(x);

x = (x + a).*b;

alpha = 0.4;
beta = -0.13;

M = 6;
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


figure(2)
plot(des,orthRep);
hold on;
plot(x,f);
hold off;
%plot(c)
