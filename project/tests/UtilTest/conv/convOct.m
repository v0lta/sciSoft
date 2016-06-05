
p1 = load('tmpP1.dat');
p2 = load('tmpP2.dat');

format short
res = conv(p1,p2);
printf("%.5f\n", res);
