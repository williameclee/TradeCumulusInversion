% Horizonatal grids
Y = linspace(-1, 1, Jt).'; % Potential latitude (S) in the paper
DY = Y(2) - Y(1);
phi = rad2deg(asin(Y));

% Vertical grids
Z = linspace(0, 1, K);
DZ = Z(2) - Z(1);
th = linspace(thetab, thetat, K).';
Dth = th(2) - th(1);
th_pad = [thetab - Dth; th; thetat + Dth];

[YY,Th_pad] = ndgrid(Y,th_pad);
[phii, ~] = ndgrid(phi, th_pad);

fac1 = eps * DY * Y;
fac2 = 2 * DY * DZ ^ 2;
f1 = fac1 .* (1 - Y .^ 2);