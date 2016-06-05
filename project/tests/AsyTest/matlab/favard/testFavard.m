clear all;

%% mon Basis
alpha = 0;
beta  = 0;

M = 3;

disp('monomial Basis')
orthPol = favard(M,alpha,beta);
omega = @(x) (1 - x) .^alpha .* (1 + x) .^beta;

%go through the rows.
for i = 1:size(orthPol,1)
    for j = 1:size(orthPol,1)
        if (i ~= j)
            %check orthogonality with respect to the other matrix entries.
            testFun = @(x) omega(x) .* polyval(orthPol(i,:),x) .* polyval(orthPol(j,:),x);
            res = integral(testFun,-1,1);
            if (res > 10^5*eps)
                disp('orth Fail')
                res
            end
            
        end
    end
    %Check the length.
    testFun = @(x) omega(x) .* polyval(orthPol(i,:),x).^2;
    res = integral(testFun,-1,1);
    if (abs(res - 1) > 10^5*eps)
        disp('fail norm')
        abs(res - 1)
    end
end


