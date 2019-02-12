function [theta, J] = LinearReg(X, y, alpha, iterTotal, Tor)
    m = length(y);
    X = [ones(m, 1), X];
    [m, n] = size(X);
    sigma = std(X);
    mu = mean(X);
    for i = 2:n
        X(:,i) = (X(:,i)-mu(i))./sigma(i);
    end 
    theta = zeros(n,1);
    J(1) = (X*theta-y)'*(X*theta-y)/(2*m);
    theta = theta-alpha*X'*(X*theta-y)/m;
    for i = 2:iterTotal
        J(i) = (X*theta-y)'*(X*theta-y)/(2*m);
        theta = theta-alpha*X'*(X*theta-y)/m;
        if(abs(J(i)-J(i-1))<Tor)
            break;
        end
    end
end
