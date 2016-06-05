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

alpha = 0.4;
beta = -0.13;

%create a polynomials;
M = 6;
lagPoly = lagrange(M);


%orthogonalize them:
orthPoly = gramSchmidt(0,0,lagPoly);
figure(1)
for i = 1:M
    plot(x,polyval(orthPoly(i,:),x));
    hold on;    
end 
hold off;

%find the p matrix:
for i = 1:M
    P(:,i) = polyval(orthPoly(i,:),x);
    plot(x,P(:,i)); hold on;
end
hold off;

c = pinv(P) * f';
orthRep = P*c;

%back transformation
x = x/b - a;

sol = load('solF.test');
sol = sort(sol);

figure(2)
plot(x,orthRep,'*');
hold on;
plot(x,f);
hold off;
%plot(c)
