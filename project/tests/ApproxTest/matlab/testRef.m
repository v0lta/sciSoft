
%define a coarse grid.
x = load('x.in');
%compute the sin on it.
f = load('fun.in');

[x,I] = sort(x);
f = f(I);

%plot
plot(x,f,'*');
hold on;


%define a finer grid.
des = load('des.in');
des = sort(des);
%des = -5:0.01:7;

%intepolate!!!
% 1. compute the weights:

lambda = zeros(length(f),1);


for i =1:length(x)
  tmp = 1;  
  for j =1:length(x)
    if (i ~= j) 
       tmp = tmp * (x(i) - x(j));       
    end
  end
  lambda(i) = 1/tmp;  
end

interpValues = zeros(length(des),1);
l = zeros(length(f),1);

for k= 1:length(des)
   for i = 1:length(f)
      tmp = 0; 
      for j = 1:length(f)
          %compute the sum of different mu's.
          tmp = tmp + lambda(j)/(des(k) - x(j));
      end
              %mu_i / sum mu_i
      l(i) = (lambda(i)/(des(k) - x(i)))/tmp;
   end
   interpValues(k) =  f*l;
end

plot(des,interpValues,'*')

figure
solF = load('solF.test');
solY= load('solY.test');
[solY,I] = sort(solY);
solF = solF(I);

plot(solY,solF,'*')

