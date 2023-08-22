%% Pressure
tbdif = thetat - thetab;
tidif = theta2 - theta1;

sint = sin(pi * (2 .* Th_pad - theta2 - theta1) / tidif);
tfac = (tidif * (Th_pad - thetab)) / (2 * tbdif);

phi_I = (phii > -phiw) & (phii < phiw);
fphi = ones(size(phii));
fphi(phi_I) = 0.5 * (1 - cos(pi * phii(phi_I) / phiw));

th_T = Th_pad > theta2;
th_C = Th_pad >= theta1 & Th_pad <= theta2;
gthet = -tfac;
gthet(th_T) = gthet(th_T) + 0.5 * tidif;
gthet(th_C) = gthet(th_C) + 0.5 * ((Th_pad(th_C) - theta1) + (tidif / 2 / pi) * sint(th_C));

Pres = (pb - sigma0 * (Th_pad - thetab) - sighat * fphi .* gthet) / pb;
Gamma = Pres .^ (kappa - 1);

pref = zeros([1, Kt]);

for k = 1:Kt
    pref(k) = pmean(Pres(:, k), Y);
end

%% M and s
X = zeros([Jt, Kt]);
X(2:2:end - 1, 2) = cp * thetab / cs;

for j = 2:2:Jt - 1
    Mp = Pres(j, 2:end - 1) .^ kappa / (kappa * alpha);
    X(j, 2:end - 1) = intdde(Mp, X(j, 2), DZ);
end

X(2:2:end - 1, 1) = X(2:2:end - 1, 3) - 2 * DZ * X(2:2:end - 1, 2) / beta;
X(2:2:end - 1, end) = X(2:2:end - 1, end - 2) + 2 * DZ * (1 - alpha) ^ kappa / (kappa * alpha);

X(1:2:end, :) = repmat(Y(1:2:end), [1, Kt]);

%% Sigma
cost = cos(pi * (2 .* Th_pad - theta2 - theta1) / tidif);
fthet = ones(size(Th_pad)) * (-tidif / (2 * tbdif));
fthet(th_C) = fthet(th_C) + 0.5 * (1 + cost(th_C));
Sigma = (sigma0 + sighat * fphi .* fthet) / sigma0;

% updateBoundary

Pres0 = Pres;
X0 = X;

clear cost fphi fthet gthet phi_I tbdif tfac th_B th_C th_T thetat tidif sint
