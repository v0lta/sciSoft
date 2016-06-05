clear all;
%% mon Basis
alpha = 0.2;
beta  = -0.2;

M = 4;
fail = 0;

combine2

omega = @(x) (1 - x) .^alpha .* (1 + x) .^beta ;

%go through the rows.
for i = 1:size(orthPol,1)
    for j = 1:size(orthPol,1)
        if (i ~= j)
            %check orthogonality with respect to the other matrix entries.
            testFun = @(x) omega(x) .* polyval(orthPol(i,:),x) .* polyval(orthPol(j,:),x);
            res = quad(testFun,-1,1);
            if (res > 10^-5)
                disp('orth Fail')
                abs(res)
                fail = 1;
            end

        end
    end
    %Check the length.
    testFun = @(x) omega(x) .* polyval(orthPol(i,:),x).^2;
    res = quad(testFun,-1,1);
    if (abs(res - 1) > 10^-5)
        disp('fail norm')
        abs(res)
        fail = 1;
    end
end

if (fail == 0)
  disp('ok')
end
