function [X, Pres, Gamma, ape, dif] = updateBoundary(X, Pres, Gamma, pref, eps, DZ, alpha, beta, kappa, Su, F1)
    X = updateBoundary_M(X, eps, DZ, Su, alpha, beta, kappa);
    X = updateBoundary_Sl(X, Su, F1);
    [Pres, Gamma, ape, dif] = correcction_Gamma(Pres, Gamma, X, pref, DZ, alpha, kappa);
    X = correcction_M(X, Pres, DZ, alpha, kappa);
    X = updateBoundary_M(X, eps, DZ, Su, alpha, beta, kappa);
end

function X = updateBoundary_M(X, eps, DZ, Su, alpha, beta, kappa)
    M_int = 2:2:size(X, 1) - 1;
    % top
    X(M_int, end) = X(M_int, end - 2) ...
        + 2 * DZ * (1 - alpha) ^ kappa / (kappa * alpha);
    % bottom
    t = (Su(1:2:end) .^ 2 - X(1:2:end, 2) .^ 2) ./ (1 - X(1:2:end, 2) .^ 2);
    t(1) = 0;
    t(end) = 0;
    t3 = eps * DZ * (t(1:end - 1) + t(2:end)) / (8 * beta);
    X(M_int, 1) = X(M_int, 1) ...
        - 2 * DZ * X(M_int, 2) / beta ...
        + t3;
end

function X = updateBoundary_Sl(X, Su, F1)
    X = updateBoundary_Sl_Iter(X, Su, F1, 1);
    X = updateBoundary_Sl_Iter(X, Su, F1, size(X, 2));
end

function X = updateBoundary_Sl_Iter(X, Su, F1, k)
    Sl_int = 3:2:size(X, 1) - 2;
    Xold = X(Sl_int, k);

    for it = 1:4
        Gkp = F1(Sl_int) .* (Xold .^ 2 - Su(Sl_int) .^ 2) ./ (1 - Xold .^ 2) ...
            + X(Sl_int + 1, k) - X(Sl_int - 1, k);
        Gkpp = F1(Sl_int) * 2 .* Xold .* (1 - Su(Sl_int) .^ 2) ./ ((1 - Xold .^ 2) .^ 2);
        X(Sl_int, k) = Xold - Gkp ./ Gkpp;
    end

    X(Sl_int, k) = Xold;

end

function [Pres, Gamma, ape, dif] = correcction_Gamma(Pres, Gamma, X, pref, DZ, alpha, kappa)
    M_int = 2:2:size(X, 1) - 1;
    Pres(M_int, 3:end - 2) = (alpha * kappa * (X(M_int, 4:end - 1) - X(M_int, 2:end - 3)) / (2 * DZ)) .^ (1 / kappa);
    ape = sum(Pres(M_int, :) .* (X(M_int + 1, :) - X(M_int - 1, :)) / 2, 1);
    dif = ape - pref;
    Pres(M_int, :) = Pres(M_int, :) - dif;
    Gamma(M_int, :) = Pres(M_int, :) .^ (kappa - 1);
    ape = sum(Pres(M_int, :) .* (X(M_int + 1, :) - X(M_int - 1, :)) / 2, 1);
    dif = ape - pref;
end

function X = correcction_M(X, Pres, DZ, alpha, kappa)
    Jt = size(X, 1);

    for j = 2:2:Jt - 1
        Mp = Pres(j, 2:end - 1) .^ kappa / (kappa * alpha);
        X(j, 2:end - 1) = intdde(Mp, X(j, 1), DZ);
    end

end
