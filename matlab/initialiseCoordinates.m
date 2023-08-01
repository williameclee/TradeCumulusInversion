% Horizonatal grids
Su = linspace(-1, 1, Jt);
DSu = Su(2) - Su(1);
Phiu = rad2deg(asin(Su));

% Vertical grids
Z = linspace(0, 1, K);
DZ = Z(2) - Z(1);
theta = linspace(thetab, thetat, K).';
theta_pad = [thetab; theta; thetat];