
actualValues = sort(load('solF.test'));
des = sort(load('des.in'));

gsinterpol    = load('gsinterpol.dat');
recurinterpol = load('recurinterpol.dat');
asyinterpol   = load('asyinterpol.dat');
gsproj     = load('gsproj.dat');
recurproj  = load('recurproj.dat');
asyproj    = load('asyproj.dat');



figure(1)
plot(des,actualValues,'*');
hold on;
plot(des,gsinterpol,'*');
plot(des,recurinterpol,'*');
plot(des,recurproj,'*')
plot(des,asyinterpol,'*');
plot(des,asyproj,'*');

legend('sampled function','gs interpol','recur interpol', 'recur proj', ...
       'asy interpol', 'asy proj'); 
