
actualValues = fliplr(load('abiclose.in'));

gsinterpol    = load('gsinterpol.dat');
recurinterpol = load('recurinterpol.dat');
asyinterpol   = load('asyinterpol.dat');
gsproj     = load('gsproj.dat');
recurproj  = load('recurproj.dat');
asyproj    = load('asyproj.dat');


x = linspace(-1,1,length(actualValues));
des = [linspace(-1,1,length(actualValues)*2) 1.1 1.2 1.3];

figure(1)
plot(x,actualValues);
hold on;
%plot(des,gsinterpol);
%plot(des,gsproj);
plot(des,recurinterpol);
plot(des(3:(end-2)),recurproj(3:(end-2)))
%plot(des,asyinterpol);
%plot(des,asyproj);

legend('stock', 'recur interpol','recur proj')

week1 = datetime('2015-01-01'):caldays(1):datetime('2015-01-02');
week2 = (week1(end)+3);
week2 = datetime('2015-01-05'):caldays(1):(week2(1)+4);
week3 = (week2(end)+3);
week3 = (week2(end)+3):caldays(1):(week3(1)+4);
week4 = (week3(end)+3);
week4 = (week3(end)+3):caldays(1):(week4(1)+4);
week5 = (week4(end)+3);
week5 = (week4(end)+3):caldays(1):(week5(1)+4);
week6 = (week5(end)+3);
week6 = (week5(end)+3):caldays(1):(week6(1)+4);
week7 = (week6(end)+3);
week7 = (week6(end)+3):caldays(1):(week7(1)+4);
week8 = (week7(end)+3);
week8 = (week7(end)+3):caldays(1):(week8(1)+4);
week9 = (week8(end)+3);
week9 = (week8(end)+3):caldays(1):(week9(1)+4);
week10 = (week9(end)+3);
week10 = (week9(end)+3):caldays(1):(week10(1)+4);
week11 = (week10(end)+3);
week11 = (week10(end)+3):caldays(1):(week11(1)+4);
week12 = (week11(end)+3);
week12 = (week11(end)+3):caldays(1):(week12(1)+4);
week13 = (week12(end)+3);
week13 = (week12(end)+3):caldays(1):(week13(1)+4);
week14 = (week13(end)+3);
week14 = (week13(end)+3):caldays(1):(week14(1)+4);
week15 = (week14(end)+3);
week15 = (week14(end)+3):caldays(1):(week15(1)+4);
time = [week1 week2 week3 week4 week5 week6 week7 week8 week9 week10 ...
        week11 week12 week13 week14 week15];


destime = time(1):hours(12):time(end);
    
    
%figure(2)
%plot(time,actualValues);
%hold on;
%plot(des,gsinterpol);
%plot(des,recurinterpol);
%plot(des,asyinterpol);
%plot(des,gsproj);
%plot(des(3:(end-2)),recurproj(3:(end-2)))
%plot(des,

