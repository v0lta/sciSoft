
data = load('ref1b.in');
format long

alpha = 0.4;
beta = -0.13;
N = length(data);

ab = r_jacobi(N,alpha,beta);
xw = gauss(N,ab);

sum = 0;
for i=1:N
    sum = sum + xw(i,2)*data(i);
end
disp (sum);
