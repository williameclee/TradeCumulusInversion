% Experiment dimension
K = 61;
J = 81;
Jt = 2 * J + 1;
Kt = K + 2;

% Physical parameters
R = 287;
cp = 1004;
kappa = R / cp;
gr = 9.81;
Ae = 6371e3;
Omega = 7.292e-5;
rom = Ae * Omega;

% Experiment parameters
thetat = 320; % theta (kelvin) at model top
thetab = 300; % theta (kelvin) at model bottom
sigma0 = 45; % sigma0 (hPa/K)

pb = 1000; % pressure (hPa) at model bottom
Dp = 58; % pressure (hPa) change across inversion
theta10 = 303; % mean theta (kelvin) in inversion
Dtheta = 4; % dtheta (kelvin) of inversion
thetac = 3; % thetac (kelvin) of inversion

phis = -30; % southern extent of inversion (degrees)
phin = 30; % northern extent of inversion (degrees)
phic = 0; % latitude (degrees) of inversion center

Dp = Dp * 100; % pressure (Pa) change across inversion
pb = pb * 100; % pressure (Pa) at model bottom
sigma0 = sigma0 * 100; % sigma0 (Pa/K)

pt = pb - sigma0 * (thetat - thetab) + Dp;
alpha = (pb - pt) / pb;
beta = thetab / (thetat - thetab);
cs = alpha * R * (thetat - thetab);
eps = 4 * rom ^ 2 / cs;

%% Create coordinates
initialiseCoordinates
% Numerical parameters
F1 = eps * DSu * Su;
F2 = 2 * DSu * DZ ^ 2;
f1 = F1 .* (1 - Su .^ 2);

%% Initialise variables
initialiseVariables

%% Visualise initial field
plot_InitialField