clear all
clc

% Experiment dimensions
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

% Default model parameters
thetat = 320; % theta (kelvin) at model top
thetab = 300; % theta (kelvin) at model bottom
sigma0 = 45 * 100; % sigma0 (Pa/K)
pb = 1000 * 100; % pressure (Pa) at model bottom

initialiseParameters_exp2

alpha = (pb - pt) / pb;
beta = thetab / (thetat - thetab);
cs = alpha * R * (thetat - thetab);
eps = 4 * rom ^ 2 / cs;

initialiseCoordinates
initialiseVariables_exp2
iterateSolution
plotFinalField_exp2
