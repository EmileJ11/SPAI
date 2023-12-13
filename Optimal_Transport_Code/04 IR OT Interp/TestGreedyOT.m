%% Simple example of PCOT object
% Created by: Aaron Geldert
% Last modified: 25 Oct 2022
clear; clc; 
% close all;
pos1 = [.75 1 1];
pos2 = [1 2 1.2];
pc1 = imageSourcePointCloud([3 5 2], pos1, 2);
pc2 = imageSourcePointCloud([3 5 2], pos2, 2);

xi = vecnorm(pos2 - pos1);
OT = PCOT(pc1, pc2, xi.^2);
scale = 900;

% figure;
subplot(121);
scatter3(pc1.pos(:,1), pc1.pos(:,2), pc1.pos(:,3), scale*pc1.mass, 'r.');
hold on;
scatter3(pc2.pos(:,1), pc2.pos(:,2), pc2.pos(:,3), scale*pc2.mass, 'b.');
% axis([-10 10 -10 10 -10 10]);
title('OT')
view([45 45]);

subplot(122);
scatter3(pc1.pos(:,1), pc1.pos(:,2), pc1.pos(:,3), scale*pc1.mass, 'r.');
hold on;
scatter3(pc2.pos(:,1), pc2.pos(:,2), pc2.pos(:,3), scale*pc2.mass, 'b.');
% axis([-10 10 -10 10 -10 10]);
title('Greedy');
view([45 45]);

for k = 0:0.01:1
    pck = OT.interpPC(k);
    pcg = OT.greedyInterp(k);
    
    subplot(121);
    h1 = scatter3(pck.pos(:,1), pck.pos(:,2), pck.pos(:,3), scale*pck.mass, 'm.');
    
    subplot(122);
    h2 = scatter3(pcg.pos(:,1), pcg.pos(:,2), pcg.pos(:,3), scale*pcg.mass, 'm.');
    
    sgtitle(['k = ' num2str(k)]);
    
    pause(0.1);
    delete(h1);
    delete(h2);
end

% it works!

