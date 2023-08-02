U = zeros(size(X));
D = zeros(size(X));
L = zeros(size(X));
G = zeros(size(X));

Sl_int = 3:2:size(X, 1) - 2;
M_int = 2:2:size(X, 1) - 1;

for n = 1:100

    for kz = 2:3

        for k = kz:2:K - 1
            G(Sl_int, k) =- (F1(Sl_int) .* ((Su(Sl_int) .^ 2 - X(Sl_int, k) .^ 2) ./ (1 - X(Sl_int, k) .^ 2))) + X(Sl_int + 1, k) - X(Sl_int - 1, k);
            D(Sl_int, k) = F1(Sl_int) * 2 .* Su(Sl_int) .* ((1 - Su(Sl_int) .^ 2) ./ (1 - X(Sl_int, k) .^ 2) .^ 2);
        end

        X(Sl_int, kz:2:K - 1) = X(Sl_int, kz:2:K - 1) - (G(Sl_int, kz:2:K - 1) ./ D(Sl_int, kz:2:K - 1));

        X(M_int, kz:2:K - 1) = X(M_int, kz:2:K - 1) - (G(M_int, kz:2:K - 1) ./ D(M_int, kz:2:K - 1));

    end
    [X, Pres, Gamma, ape, dif] = updateBoundary(X, Pres, Gamma, pref, eps, DZ, alpha, beta, kappa, Su, F1);

end
