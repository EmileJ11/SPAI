%% SDM implementation
% Created by: Aaron Geldert and Nils Meyer-Kahlen
% Last modified: 25 May 2022

clear; clc; close all;

% Load ambisonic SRIR data
filename = 'm5_ls1.sofa';
D = SOFAload(filename);
fs = D.Data.SamplingRate; 
c = 343;

threshold = 32;

% Convert SRIRs to full SDM point clouds
srir1 = squeeze(D.Data.IR(15, :, 301:4000))';
srir2 = squeeze(D.Data.IR(17, :, 301:4000))';
SDM1 = processSDM(srir1, fs, 1, threshold);
SDM2 = processSDM(srir2, fs, 1, threshold);

% Visualize IR and point cloud
scale = 900;

%% Interpolation
xi = 0.1;
k = 0.5;
sRatio = 0.8;
SDMint = interpolateSDM(SDM1, SDM2, k, xi, sRatio);

% Render to loudspeakers
tdsgn = getTdesign(6);
[doas1, p1] = fillPointCloud(SDM1);
[doasInt, pInt] = fillPointCloud(SDMint);
[doas3, p3] = fillPointCloud(SDM2);
x1 = encodeToLs(doas1, tdsgn, p1);
x2Int = encodeToLs(doasInt, tdsgn, pInt);
x3 = encodeToLs(doas3, tdsgn, p3);

figure; 
subplot(311); plot(x1 + (0:23)); xlim([0 1e3]);
subplot(312); plot(x2Int + (0:23)); xlim([0 1e3]);
subplot(313); plot(x3 + (0:23)); xlim([0 1e3]);

% NOTES:
% Interpolation on cartesian coordinates means that things can pass by the
% origin - better to interpolate spherically? i.e. correct the distance?


%% Visualize interpolation
axlim = [-5 5 -10 5 -5 5];
scale = 2;

figure;
plotSDM(subplot(311), SDM1, '1', axlim, scale);
plotSDM(subplot(312), SDMint, '2 interp', axlim, scale);
plotSDM(subplot(313), SDM2, '3', axlim, scale);



