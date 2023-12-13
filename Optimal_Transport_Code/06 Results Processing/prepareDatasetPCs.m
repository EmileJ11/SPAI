%% Prepare the test datasets
% Created by: Aaron Geldert
% Last modified: 10 Oct 2022

% 1. Define the geometry, source, reciever positions
% 2. Define the ISM model parameters
% 3. Process the ISM model to generate PCs for all positions
% 4. Save all data

tic;
%% Cuboid, MR, SS, 1
clear; close all;

pathout = '/Users/aarongeldert/Documents/MATLAB/THESIS/ISM_Data/';

params.filename     = 'cuboid.txt';
params.srcPos       = [6 3.6 1.6];
params.rcvPos       = [2 4 1.8];
params.movDir       = [1.5 -1 0];
params.movInterval  = 0.10;
params.numMovPos    = 21;
params.moving       = 'rcv';
params.maxIsOrder   = 1;
params.fs           = 48e3;
params.doPlot       = 0;

if ~exist('globalSurfaceCount','var')
    global globalSurfaceCount
end
globalSurfaceCount = 1;
[PCs, config] = runISM(params);
fileout = [pathout 'Cuboid_moving_' params.moving '_IS_' num2str(params.maxIsOrder) '.mat'];
save(fileout, 'PCs', 'config', 'params');

%% Cuboid, MR, SS, 5
clear; close all;
pathout = '/Users/aarongeldert/Documents/MATLAB/THESIS/ISM_Data/';

params.filename     = 'cuboid.txt';
params.srcPos       = [6 3.6 1.6];
params.rcvPos       = [2 4 1.8];
params.movDir       = [1.5 -1 0];
params.movInterval  = 0.10;
params.numMovPos    = 21;
params.moving       = 'rcv';
params.maxIsOrder   = 5;
params.fs           = 48e3;
params.doPlot       = 0;

if ~exist('globalSurfaceCount','var')
    global globalSurfaceCount
end
globalSurfaceCount = 1;
[PCs, config] = runISM(params);
fileout = [pathout 'Cuboid_moving_' params.moving '_IS_' num2str(params.maxIsOrder) '.mat'];
save(fileout, 'PCs', 'config', 'params');

%% Cuboid, SR, MS, 1
clear; close all;

pathout = '/Users/aarongeldert/Documents/MATLAB/THESIS/ISM_Data/';

params.filename     = 'cuboid.txt';
params.srcPos       = [6 3.6 1.6];
params.rcvPos       = [2 4 1.8];
params.movDir       = [-4 -1 0];
params.movInterval  = 0.10;             % different!
params.numMovPos    = 21;
params.moving       = 'src';
params.maxIsOrder   = 1;
params.fs           = 48e3;
params.doPlot       = 0;

if ~exist('globalSurfaceCount','var')
    global globalSurfaceCount
end
globalSurfaceCount = 1;
[PCs, config] = runISM(params);
fileout = [pathout 'Cuboid_moving_' params.moving '_IS_' num2str(params.maxIsOrder) '.mat'];
save(fileout, 'PCs', 'config', 'params');

%% Cuboid, SR, MS, 5
clear; close all;
pathout = '/Users/aarongeldert/Documents/MATLAB/THESIS/ISM_Data/';

params.filename     = 'cuboid.txt';
params.srcPos       = [6 3.6 1.6];
params.rcvPos       = [2 4 1.8];
params.movDir       = [-4 -1 0];
params.movInterval  = 0.10;
params.numMovPos    = 21;
params.moving       = 'src';
params.maxIsOrder   = 5;
params.fs           = 48e3;
params.doPlot       = 0;

if ~exist('globalSurfaceCount','var')
    global globalSurfaceCount
end
globalSurfaceCount = 1;
[PCs, config] = runISM(params);
fileout = [pathout 'Cuboid_moving_' params.moving '_IS_' num2str(params.maxIsOrder) '.mat'];
save(fileout, 'PCs', 'config', 'params');

%% General, MR, SS, 1
clear; close all;

pathout = '/Users/aarongeldert/Documents/MATLAB/THESIS/ISM_Data/';

params.filename     = 'lofted.txt';
params.srcPos       = [1.95 5.2 1.85];
params.rcvPos       = [3.9 1.4 1.6];
params.movDir       = [-2 1 0];
params.movInterval  = 0.10;
params.numMovPos    = 21;
params.moving       = 'rcv';
params.maxIsOrder   = 1;
params.fs           = 48e3;
params.doPlot       = 0;

if ~exist('globalSurfaceCount','var')
    global globalSurfaceCount
end
globalSurfaceCount = 1;
[PCs, config] = runISM(params);
fileout = [pathout 'General_moving_' params.moving '_IS_' num2str(params.maxIsOrder) '.mat'];
save(fileout, 'PCs', 'config', 'params');

%% General, MR, SS, 5
clear; close all;
pathout = '/Users/aarongeldert/Documents/MATLAB/THESIS/ISM_Data/';

params.filename     = 'lofted.txt';
params.srcPos       = [1.95 5.2 1.85];
params.rcvPos       = [3.9 1.4 1.6];
params.movDir       = [-2 1 0];
params.movInterval  = 0.10;
params.numMovPos    = 21;
params.moving       = 'rcv';
params.maxIsOrder   = 5;
params.fs           = 48e3;
params.doPlot       = 1;

if ~exist('globalSurfaceCount','var')
    global globalSurfaceCount
end
globalSurfaceCount = 1;
[PCs, config] = runISM(params);
% fileout = [pathout 'General_moving_' params.moving '_IS_' num2str(params.maxIsOrder) '.mat'];
% save(fileout, 'PCs', 'config', 'params');

%% General, SR, MS, 1
clear; close all;

pathout = '/Users/aarongeldert/Documents/MATLAB/THESIS/ISM_Data/';

params.filename     = 'lofted.txt';
params.srcPos       = [1.95 5.2 1.85];
params.rcvPos       = [3.9 1.4 1.6];
params.movDir       = [-1 -4 0];
params.movInterval  = 0.10;             % different!
params.numMovPos    = 21;
params.moving       = 'src';
params.maxIsOrder   = 1;
params.fs           = 48e3;
params.doPlot       = 0;

if ~exist('globalSurfaceCount','var')
    global globalSurfaceCount
end
globalSurfaceCount = 1;
[PCs, config] = runISM(params);
fileout = [pathout 'General_moving_' params.moving '_IS_' num2str(params.maxIsOrder) '.mat'];
save(fileout, 'PCs', 'config', 'params');

%% General, SR, MS, 5
clear; close all;
pathout = '/Users/aarongeldert/Documents/MATLAB/THESIS/ISM_Data/';

params.filename     = 'lofted.txt';
params.srcPos       = [1.95 5.2 1.85];
params.rcvPos       = [3.9 1.4 1.6];
params.movDir       = [-2 1 0];
params.movInterval  = 0.10;
params.numMovPos    = 21;
params.moving       = 'src';
params.maxIsOrder   = 5;
params.fs           = 48e3;
params.doPlot       = 0;

if ~exist('globalSurfaceCount','var')
    global globalSurfaceCount
end
globalSurfaceCount = 1;
[PCs, config] = runISM(params);
fileout = [pathout 'General_moving_' params.moving '_IS_' num2str(params.maxIsOrder) '.mat'];
save(fileout, 'PCs', 'config', 'params');

toc;