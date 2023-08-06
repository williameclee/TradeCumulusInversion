function Y = fexp(X)
    gamma = 1/2 * exp(2) * log(2);
    Y = exp(-gamma * exp(1 ./ (X - 1)) ./ X);
end
