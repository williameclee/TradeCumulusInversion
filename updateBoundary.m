X = updateBoundary_M(X, Y, DZ, alpha, beta, eps, kappa);
X = updateBoundary_s(X, Y, fac1);
[Pres, Gamma, pmean, dif] = correction_Gamma(X, DZ, Pres, Gamma, pref, alpha, kappa);
X = correction_M(X, DZ, Pres, alpha, kappa);
X = updateBoundary_M(X, Y, DZ, alpha, beta, eps, kappa);

function X = updateBoundary_M(X, Y, DZ, alpha, beta, eps, kappa)
    M_int = 2:2:size(X, 1) - 1;
    % top
    X(M_int, end) = X(M_int, end - 2) ...
        + 2 * DZ * (1 - alpha) ^ kappa / (kappa * alpha);
    % bottom
    t = (Y(1:2:end) .^ 2 - X(1:2:end, 2) .^ 2) ./ (1 - X(1:2:end, 2) .^ 2);
    t(1) = 0;
    t(end) = 0;
    t3 = eps * DZ * (t(1:end - 1) + t(2:end)) / (8 * beta);
    X(M_int, 1) = X(M_int, 3) ...
        - 2 * DZ * X(M_int, 2) / beta ...
        + t3;
end

function X = updateBoundary_s(X, Y, fac1)
    X = updateBoundary_s_Iter(X, Y, fac1, 1);
    X = updateBoundary_s_Iter(X, Y, fac1, size(X, 2));
end

function X = updateBoundary_s_Iter(X, Y, fac1, k)
    Sl_int = 3:2:size(X, 1) - 2;

    if k == 1
        Xold = X(Sl_int, k + 1);
    else
        Xold = X(Sl_int, k - 1);
    end

    for it = 1:4
        Gkp = fac1(Sl_int) ...
            .* (Xold .^ 2 - Y(Sl_int) .^ 2) ./ (1 - Xold .^ 2) ...
            + X(Sl_int + 1, k) - X(Sl_int - 1, k);
        Gkpp = fac1(Sl_int) * 2 .* Xold ...
            .* (1 - Y(Sl_int) .^ 2) ./ ((1 - Xold .^ 2) .^ 2);
        X(Sl_int, k) = Xold - Gkp ./ Gkpp;
        Xold = X(Sl_int, k);
    end
end

function [Pres, Gamma, pmean, dif] = correction_Gamma(X, DZ, Pres, Gamma, pref, alpha, kappa)
    M_int = 2:2:size(X, 1) - 1;

    Pres(M_int, 3:end - 2) = (alpha * kappa * (X(M_int, 4:end - 1) - X(M_int, 2:end - 3)) / (2 * DZ)) .^ (1 / kappa);
    pmean = sum(Pres(M_int, :) .* (X(M_int + 1, :) - X(M_int - 1, :)), 1) / 2;
    dif = pmean - pref;
    Pres(M_int, 2:end - 1) = Pres(M_int, 2:end - 1) - dif(2:end - 1);
    Gamma(M_int, 2:end - 1) = Pres(M_int, 2:end - 1) .^ (kappa - 1);
    pmean = sum(Pres(M_int, :) .* (X(M_int + 1, :) - X(M_int - 1, :)), 1) / 2;
    dif = pmean - pref;
end

function X = correction_M(X, DZ, Pres, alpha, kappa)
    Jt = size(X, 1);

    for j = 2:2:Jt - 1
        Mp = Pres(j, 2:end - 1) .^ kappa / (kappa * alpha);
        X(j, 2:end - 1) = intdde(Mp, X(j, 2), DZ);
    end

end
