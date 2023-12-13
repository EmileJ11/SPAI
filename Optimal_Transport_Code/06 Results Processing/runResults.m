% Determinations:
clc; close all; clear;
% maxIsOrder: 3 fast, 4 complete
% good distance: 1m
% good noiseLvl: 0.05;
% good mapTol: 1.05 (5% tolerance on the movement distance)

% params.order = 3;
% params.dist = 0.5; 
% params.noise = 0.05; % not in use
% % params.paths = 200;
% params.paths = 100;
% params.nTarget = 16; % not in use 
% % statisticalCaseResults(params.order, params.dist, params.noise, params.paths, params.nTarget, 1);
% % clearvars -except params; close all;
% % statisticalCaseResults(params.order, params.dist, params.noise, params.paths, params.nTarget, 1.2);
% % clearvars -except params; close all;
% statisticalCaseResults(params.order, params.dist, params.noise, params.paths, params.nTarget, 1.01, 1, 0);
% % clearvars -except params; close all;
% 

%%
% statisticalCaseResults(params.order, params.dist, params.paths, params.mapTol, 1, 0);

clc; close all; clear;

dists = [2 1 0.5 4 0.25 0.125];
tic;
for distInd = 1:numel(dists)
    statisticalCaseResults('cuboid.txt', 3, dists(distInd), 50, 1.05, 1, 0);
    statisticalCaseResults('canted.txt', 3, dists(distInd), 50, 1.05, 1, 0);
    statisticalCaseResults('trapezoidal.txt', 3, dists(distInd), 50, 1.05, 1, 0);
end
toc;


