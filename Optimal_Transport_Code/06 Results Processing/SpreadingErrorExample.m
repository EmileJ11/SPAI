%% one example of error


LD = load('/Users/aarongeldert/Documents/MATLAB/THESIS/ISM_Data/Cuboid_moving_rcv_IS_1.mat');

pc1 = LD.PCs(5);
pc2 = LD.PCs(13);
pcTrue = LD.PCs(9);
mapTol = 2;
xi = mapTol * 8 * LD.config.movInterval;
OT = PCOT(pc1, pc2, xi);

k = 0.5;
[PCot, ~, ~] = OT.interpPC(k);
PClin = combinePCs(pc1, pc2, k);

fs = 48e3;
doPlot = 1;
Ep = getIRerror(PCot, pcTrue, fs, doPlot);
El = getIRerror(PClin, pcTrue, fs, doPlot);