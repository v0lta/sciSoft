
%define a coarse grid.
x = -2*pi:0.5:2*pi;
%compute the sin on it.
f = cos(x);

%plot
plot(x,f);
hold on;


%define a finer grid.
des = -pi:0.1:pi;

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

%plot(des,interpValues)

load('test.out')
plot(des,test)
