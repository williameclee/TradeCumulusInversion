% Experiment parameters
sighat = -44 * 100; % inversion strength (Pa/K)
theta1 = 302; % theta_1 (K) inversion bottom
theta2 = 307; % theta_2 (K) inversion top

phis = -30; % southern extent of plot (degrees)
phin = 30; % northern extent of plot (degrees)
phic = 0; % latitude (degrees) of inversion centre
phiw = 20; % width of region between inversions (degrees)

pt = pb - sigma0 * (thetat - thetab);
