clear all;
M = 6
polyMat = [  0.00  0.00  0.00  0.00  0.00  0.64
             0.00  0.00  0.00  0.00  1.32  0.55
             0.00  0.00  0.00  2.67  1.25 -0.57
             0.00  0.00  5.37  2.61 -2.41 -0.57
             0.00 10.77  5.34 -7.44 -2.44  0.56
             10.48  5.27 -9.81 -3.68  1.69  0.28 ];
x = [1, -1];
for i = 1:M
    P(:,i) = polyval(polyMat(i,:),x);
end
P