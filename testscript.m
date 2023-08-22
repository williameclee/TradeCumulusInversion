D = zeros([size(X(1:2:end,:),1),1]);
for j = 1:2:Jt
    j_d = (j+1)/2;
    D(j_d) = fac1(j) * (2 * X0(j, k)) * (1 - Y(j) ^ 2) / (1 - X0(j, k) ^ 2) ^ 2;
end
plot(D)