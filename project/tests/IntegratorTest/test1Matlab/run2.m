


x = -1:0.1:1;
data = x.^2;

alpha = 0;
beta =  0;
N = length(data);

ab = r_jacobi(N,alpha,beta);
xw = gauss(N,ab);

sum = 0;
for i=1:N
    sum = sum + xw(i,2)*data(i);
end
sum
trapz = trapz(x,data)