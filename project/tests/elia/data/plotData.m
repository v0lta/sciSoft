
actualValues = load('elia1115.in');
x = linspace(-1,1,length(actualValues));

des = linspace(-1,1,length(actualValues)*2);

gsinterpol    = load('gsinterpol.dat');
recurinterpol = load('recurinterpol.dat');
asyinterpol   = load('asyinterpol.dat');
gsproj     = load('gsproj.dat');
recurproj  = load('recurproj.dat');
asyproj    = load('asyproj.dat');



figure(1)
plot(x,actualValues,'*');
hold on;
plot(des,gsinterpol,'*');
plot(des,recurinterpol,'*');
plot(des,recurproj,'*')
%plot(des,asyinterpol,'*');
%plot(des,asyproj,'*');

legend('sampled function','gs interpol','recur interpol', 'recur proj', ...
       'asy interpol', 'asy proj'); 
