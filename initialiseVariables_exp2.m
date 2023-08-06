fac1 = eps * DY * Y;
fac2 = 2 * DY * DZ ^ 2;
f1 = fac1 .* (1 - Y .^ 2);

%% Theta
phi_N = (phi <= phin) & (phi > phic);
phi_S = (phi >= phis) & (phi < phic);
phi_C = phi == phic;
th1 = ones([Jt, 1]) * theta10;
th1(phi_N) = theta10 + Dtheta * fexp((phin - phi(phi_N)) / (phin - phic));
th1(phi_S) = theta10 + Dtheta * fexp((phi(phi_S) - phis) / (phic - phis));
th1(phi_C) = theta10 + Dtheta;
th2 = th1 + thetac;

clear phi_N phi_S phi_C

%% Pressure
[th1_mesh, ~] = ndgrid(th1, th_pad);
[th2_mesh, ~] = ndgrid(th2, th_pad);
th_T = Th_pad >= th2_mesh;
th_C = Th_pad > th1_mesh & Th_pad < th2_mesh;
Pres = pb - sigma0 * (Th_pad - thetab);
Pres(th_T) = Pres(th_T) + Dp;
Pres(th_C) = Pres(th_C) ...
    + Dp * fexp((Th_pad(th_C) - th1_mesh(th_C)) / thetac);
Pres = Pres / pb;
pref = zeros([1,Kt]);
for k = 1:Kt
    pref(k) = pmean(Pres(:, k), Y);
end
Gamma = Pres .^ (kappa - 1);

%% M and Phil
X = zeros([Jt, Kt]);
X(2:2:end - 1, 1) = cp * thetab / cs;

for j = 2:2:Jt - 1
    Mp = Pres(j, 2:end - 1) .^ kappa / (kappa * alpha);
    X(j, 2:end - 1) = intdde(Mp, X(j, 1), DZ);
end

X(2:2:end - 1, 1) = X(2:2:end - 1, 3) - 2 * DZ * X(2:2:end - 1, 2) / beta;
X(2:2:end - 1, end) = X(2:2:end - 1, end - 2) + 2 * DZ * (1 - alpha) ^ kappa / (kappa * alpha);

X(1:2:end, :) = repmat(Y(1:2:end), [1, Kt]);
% X(3:2:end - 2) = X(3:2:end - 2) + (rand(size(X(3:2:end - 2))) - 0.5) * 0.1;

%% Sigma
tinc = 0.1;
Sigma = ones([Jt, Kt]);

for j = 1:Jt

    for k = 2:Kt - 1

        if th_C(j, k) ~= 1
            continue
        end

        tp = th(k) + tinc;
        tm = th(k) - tinc;
        fp = (tp - th1(j)) / (th2(j) - th1(j));
        fm = (tm - th1(j)) / (th2(j) - th1(j));
        dfdt = (fp - fm) / (2 * tinc);
        Sigma(j, k) = (sigma0 - Dp * dfdt) / sigma0;
    end

end

clear th1_mesh th2_mesh th_T th_C

updateBoundary

Pres0 = Pres;
X0 = X;