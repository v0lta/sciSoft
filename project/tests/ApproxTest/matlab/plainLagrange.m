
%define a coarse grid.
x = -2*pi:0.5:2*pi;
%compute the sin on it.
f = exp(x);

%define a finer grid.
des = -pi/16;

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
          %compute the sum of different mus.
          tmp = tmp + lambda(j)/(des(k) - x(j));
      end
              %mu_i / sum mu_i
      l(i) = (lambda(i)/(des(k) - x(i)))/tmp;
   end
   interpValues(k) =  f*l;
end

disp(interpValues)
