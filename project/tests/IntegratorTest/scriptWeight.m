data = load('tests/ref1b.in');

alpha = 0.4;
beta = -0.13;


N = length(data)+1;
a = alpha;
b = beta;


nu=(b-a)/(a+b+2);
mu=2^(a+b+1)*gamma(a+1)*gamma(b+1)/gamma(a+b+2);

if N==1, ab=[nu mu];
    return
end

N=N-1;
n=1:N;
nab=2*n+a+b;
%size A: N + 1;
A = [nu (b^2-a^2)*ones(1,N)./(nab.*(nab+2))];
n = 2:N; nab = nab(n);
B1 =  4*(a+1)*(b+1)/((a+b+2)^2*(a+b+3));
B  =  4*(n+a).*(n+b).*n.*(n+a+b)./((nab.^2).*(nab+1).*(nab-1));
ab =  [A' [mu; B1; B']];

disp (ab)
