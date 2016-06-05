
function lagPoly = lagrange(deg)

x = linspace(-1,1,deg); 


for i =1:deg
  denom = 1;
  nom = 1;
  for j =1:length(x)
    if (i ~= j) 
       denom = denom * 1/(x(i) - x(j));
       %[1 -x(j)]
       %nom
       nom = conv(nom, [1 -x(j)]);
       %disp ('---------')
    end
  end
  lagPoly(i,:) = nom .* denom;
end


