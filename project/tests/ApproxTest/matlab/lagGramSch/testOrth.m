clear all;

%% mon Basis
alpha = 0;
beta  = 0;

M = 3;


disp('monomial Basis')
monPoly = flipud(eye(M));
orthPol = gramSchmidt( alpha,beta,monPoly );

omega = @(x) (1 - x) .^alpha .* (1 + x) .^beta ;

%go through the rows.
for i = 1:size(orthPol,1)
    for j = 1:size(orthPol,1)
        if (i ~= j)
            %check orthogonality with respect to the other matrix entries.
            testFun = @(x) omega(x) .* polyval(orthPol(i,:),x) .* polyval(orthPol(j,:),x);
            res = integral(testFun,-1,1);
            if (res > eps)
                disp('orth Fail')
                res
            end
            
        end
    end
    %Check the length.
    testFun = @(x) omega(x) .* polyval(orthPol(i,:),x).^2;
    res = integral(testFun,-1,1);
    if (abs(res - 1) > eps)
        disp('fail norm')
        abs(res - 1)
    end
end


%% Lagrange Basis
disp('Lagrange Basis')
lagPoly = lagrange(M);
orthPol = gramSchmidt( alpha,beta,lagPoly );


omega = @(x) (1 - x) .^ alpha .* (1 + x) .^beta ;

for i = 1:size(orthPol,1)
    for j = 1:size(orthPol,1)
        if (i ~= j)
            %check orthogonality with respect to the other matrix entries.
            testFun = @(x) omega(x) .* polyval(orthPol(i,:),x) .* polyval(orthPol(j,:),x);
            res = integral(testFun,-1,1);
            if (res > eps)
                disp('orth Fail')
                res
            end
            
        end
    end
    %Check the length.
    testFun = @(x) omega(x) .* polyval(orthPol(i,:),x).^2;
    res = integral(testFun,-1,1);
    if (abs(res - 1) > eps)
        disp('fail norm')
        abs(res - 1)
    end
end
