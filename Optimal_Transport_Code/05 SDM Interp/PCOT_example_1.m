%% SDM interpolation with PCOT
% Created by: Aaron Geldert
% Last modified: 31 May 2022

clear; clc; close all;

% Load ambisonic SRIR data - Nils MK
filename = 'm5_ls1.sofa';
D = SOFAload(filename);
fs = D.Data.SamplingRate; 
c = 343;

srir = permute(D.Data.IR, [3 2 1]) ; % [samples x chan x measPos]

% normalize IR
srir = srir./max(abs(srir(:)));

% Create SDMPC data structs
samps = 301:4000;
numSamples = numel(samps);

% maybe try positions between 22-28?
measPos = 1:3;
s1 = SDMPC(fs, srir(samps, :, measPos(1), :), 'foo'); % 'foo' just marks that we use B-format data
s2 = SDMPC(fs, srir(samps, :, measPos(2), :), 'foo');
s3 = SDMPC(fs, srir(samps, :, measPos(3), :), 'foo');

%% Example 2: SDM with closed walls

clear; clc; close all;

% Otto Puomio SDM dataset
filename = 'SDM_r1.sofa';
D = SOFAload(filename);
fs = D.Data.SamplingRate; 
c = 343;

rir = squeeze(D.Data.IR).'; % [samples x measPos]
doa = 0.1* D.EmitterPosition; % [samples x doa x measPos] 
% ^ CONVERSION FROM CM TO M - maybe??

% normalize IR
rir = rir./max(abs(rir(:)));

samps = 1:10000;
numSamples = numel(samps);
s1 = SDMPC(fs, squeeze(rir(samps, 3)), doa(samps, :, 1));
s2 = SDMPC(fs, squeeze(rir(samps, 4)), doa(samps, :, 2));
s3 = SDMPC(fs, squeeze(rir(samps, 4)), doa(samps, :, 3));

%% Select subset of points

% OLD METHOD: threshold & num
threshold = -60;
% maxNum = 3;
% s1.selectPoints(threshold, maxNum);
% s2.selectPoints(threshold, maxNum);
% s3.selectPoints(threshold, maxNum);

% NEW METHOD: nearly balanced energy dists
nTarget = 100;
[s1, s3] = SDMPC.selectBalancedPCs(s1, s3, nTarget);

s2.selectPoints(threshold, nTarget);

pc1 = s1.getPC();
pc2 = s2.getPC();
pc3 = s3.getPC();

df1 = s1.getDF(); % TODO check out DF thingy
df3 = s3.getDF();

% Prep the PCOT object
maxDist = 4; 
tic;
ot13 = PCOT(s1.PC, s3.PC, maxDist^2);
toc;
[pctm, pcsm] = ot13.interpPC(0.0);

figure(1);
imagesc(ot13.Tx); title('Optimal Transport Matrix Tx');

%% Loudspeaker plot
% Interpolate
kTest = 0.5;
[pctm, pcsm] = ot13.interpPC(kTest);

pck.pos = vertcat(pctm.pos, pcsm.pos);
pck.mass = vertcat(pctm.mass, pcsm.mass);
pck.n = length(pck.mass);

% Render to loudspeakers
tdsgn = getTdesign(6);

scale = 7;

x1 = scale * encodePCtoLS(pc1, tdsgn, fs);
x2 = scale * encodePCtoLS(pc2, tdsgn, fs);
xk = scale * encodePCtoLS(pck, tdsgn, fs);
x3 = scale * encodePCtoLS(pc3, tdsgn, fs);

% Render diffuse fields
% d1 = scale * encodeToLs(df1.pos, tdsgn, df1.mass);
% d3 = scale * encodeToLs(df3.pos, tdsgn, df3.mass);
% dk = (1-kTest)*d1 + kTest*d3;
d1 = 0; d3 = 0; dk = 0; % dummy solution, no diffuse field!

figure(2); 
subplot(411); plot(x1 + d1 + (0:23)); xlim([1 1e3]); title(['Pos ' num2str(measPos(1))]);
subplot(412); plot(x2 + (0:23)); xlim([1 1e3]); title(['Pos ' num2str(measPos(2)) ' meas']);
subplot(413); plot(xk + dk + (0:23)); xlim([1 1e3]); title(['Pos ' num2str(measPos(2)) ' interp: k = ' num2str(kTest)]);
subplot(414); plot(x3 + d3 + (0:23)); xlim([1 1e3]); title(['Pos ' num2str(measPos(3))]);

%% Animated
pause;
for kk = 0:0.01:1
    subplot(413);
    [pctm, pcsm] = ot13.interpPC(kk);
    pck.pos = vertcat(pctm.pos, pcsm.pos);
    pck.mass = vertcat(pctm.mass, pcsm.mass);
    pck.n = length(pck.mass);
    [doask, pk] = fillPC(pck, fs, numSamples); % debug!
    xk = scale * encodeToLs(doask, tdsgn, pk);
    dk = (1-kk)*d1 + kk*d3;
    plot(xk + dk + (0:23)); xlim([0 1e3]);
    title(['Pos ' num2str(measPos(2)) ': k = ' num2str(kk)]); 
    ylim([0 25]);
    pause(0.06);
end

%% 3D Plot with animation
axlim = [-5 5 -10 5 -5 5];
scale = 100;

figure(3);
subplot(221);
scatter3(0,0,0,'*'); % origin  (listener position)
hold on; grid on;
scatter3(pc1.pos(:,1), pc1.pos(:,2), pc1.pos(:,3),...
    pc1.mass * scale, 'ko','filled');
view([0 90]);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('PC1');
axis(axlim);

subplot(223);
scatter3(0,0,0,'*'); % origin  (listener position)
hold on; grid on;
scatter3(pc2.pos(:,1), pc2.pos(:,2), pc2.pos(:,3),...
    pc2.mass * scale, 'ko','filled');
view([0 90]);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('PC2');
axis(axlim);


subplot(224);
scatter3(0,0,0,'*'); % origin  (listener position)
hold on; grid on;
scatter3(pc3.pos(:,1), pc3.pos(:,2), pc3.pos(:,3),...
    pc3.mass * scale, 'ko','filled');
view([0 90]);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('PC3');
axis(axlim);


subplot(222);
pause;
for kk = 0:0.01:1
    [pctm, pcsm] = ot13.interpPC(kk);
    pck.pos = vertcat(pctm.pos, pcsm.pos);
    pck.mass = vertcat(pctm.mass, pcsm.mass);
    pck.n = length(pck.mass);
    
    scatter3(pck.pos(:,1), pck.pos(:,2), pck.pos(:,3),...
        eps+pck.mass * scale, 'ko','filled');
    view([0 90]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(['PC, k = ' num2str(kk)]);
    axis(axlim);
    pause(0.05);
end

%% Notes:

% Are interpolation k=0 and k=1 really matching?
% Note that point cloud scatter data will not combine several points that overlap




