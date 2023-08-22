% Experiment parameters
Dp = 58 * 100; % pressure (hPa) change across inversion
theta10 = 303; % mean theta (kelvin) in inversion
Dtheta = 4; % dtheta (kelvin) of inversion
thetac = 3; % thetac (kelvin) of inversion

phis = -30; % southern extent of inversion (degrees)
phin = 30; % northern extent of inversion (degrees)
phic = 0; % latitude (degrees) of inversion centre

pt = pb - sigma0 * (thetat - thetab) + Dp;