load("workspace_just_tracksFinal.mat")
dt    = 0.04;
Label = "DATANAME";

R     = 1/6

% import tracksFinal, somehow,
% then run this command:

[ life, dens, lumen, Diffdat_p, D_observations] = polytrack( tracksFinal, Label, dt, R);

