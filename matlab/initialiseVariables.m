%% Theta
theta1 = ones([Jt, 1]) * theta10;

phi_N = (Phiu <= phin) & (Phiu > phic);
phi_S = (Phiu >= phis) & (Phiu < phic);
phi_C = Phiu == phic;
theta1(phi_N) = theta10 + Dtheta * fexp((phin - Phiu(phi_N)) / (phin - phic));
theta1(phi_S) = theta10 + Dtheta * fexp((Phiu(phi_S) - phis) / (phic - phis));
theta1(phi_C) = theta10 + Dtheta;
theta2 = theta1 + thetac;

clear phi_N phi_S phi_C

%% Pressure
[theta1_mesh, Theta_pad] = ndgrid(theta1, theta_pad);
[theta2_mesh, ~] = ndgrid(theta2, theta_pad);
theta_T = Theta_pad >= theta2_mesh;
theta_C = Theta_pad > theta1_mesh & Theta_pad < theta2_mesh;
Pres = pb - sigma0 * (Theta_pad - thetab);

Pres(theta_T) = Pres(theta_T) + Dp;
Pres(theta_C) = Pres(theta_C) ...
    + Dp * fexp((Theta_pad(theta_C) - theta1_mesh(theta_C)) / thetac);
Pres = Pres / pb;
pref = mean(Pres, 1);
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

X(1:2:end, :) = repmat(Su(1:2:end), [1, Kt]);

%% Sigma
tinc = 0.1;
Sigma = ones([Jt, Kt]);

for j = 1:Jt

    for k = 2:Kt - 1

        if theta_C(j, k) ~= 1
            continue
        end

        tp = theta(k) + tinc;
        tm = theta(k) - tinc;
        fp = (tp - theta1(j)) / (theta2(j) - theta1(j));
        fm = (tm - theta1(j)) / (theta2(j) - theta1(j));
        dfdt = (fp - fm) / (2 * tinc);
        Sigma(j, k) = (sigma0 - Dp * dfdt) / sigma0;
    end

end

clear theta1_mesh theta2_mesh theta_T theta_C