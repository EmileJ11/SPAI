%% Make Geometrical Acoustics RIR plot animation
clear; clc; close all;
c = 343;
fs = 200*c; % 68600;

% Geometry
roomDims = [7 4 3];
srcPos = [.72 0.9 .5];

maxDist = 36;

% Determine image sources
PC = imageSourcePointCloud(roomDims, srcPos, 12, 0); % no atten calc of mass
% remove the Z direction
inds2D = PC.pos(:,3)==srcPos(3);
% remove the ones outside the fully defined distance
dists = vecnorm(PC.pos,2,2);
indsDef = dists<maxDist;

indsValid = inds2D & indsDef;

PC.pos = PC.pos(indsValid,:);
PC.mass = PC.mass(indsValid);
PC.n = numel(PC.mass);

% Receiver stuff
rcvPos = [-1.5 -0.2 .5];
PCrel.pos = PC.pos - rcvPos;
PCrel.mass = PC.mass./vecnorm(PCrel.pos, 2, 2);
PCrel.n = PC.n;
[rir, tVec] = monoIRfromPC(PCrel, maxDist/c, fs);

gifname = 'rir_idealwaves.gif';
tgif = 1/30;

% Plot setup
figure(1);
% subplot(4,1,1:3);
% scatter(PC.pos(:,1), PC.pos(:,2), 70, 'ko','filled');
scatter(PC.pos(:,1), PC.pos(:,2), 10, 'ko','filled');
hold on;
% scatter(rcvPos(1), rcvPos(2), 70, 'ko','LineWidth',1.5);
scatter(rcvPos(1), rcvPos(2), 30, 'ko','LineWidth',1.5);
% line([rcvPos(1)+0.08 rcvPos(1)+0.08], [rcvPos(2)-0.1 rcvPos(2)+0.1],'Color','black','LineWidth',1.5); 
rectangle('Position',[-0.5*roomDims(1:2), roomDims(1:2)]);
view(0,90);
grid off;
axis equal;
hr = scatter(PC.pos(1), PC.pos(2),'k.','MarkerEdgeAlpha',0);
axis([-0.5*roomDims(1) 0.5*roomDims(1) -0.5*roomDims(2) 0.5*roomDims(2)]);
xticks({}); yticks({}); 
% xlabel('X'); ylabel('Y');
title('\textbf{Ideal room impulse response, with virtual sources}');

subplot(4,1,4); cla;
hs = stem(tVec, db(rir),'k','filled','BaseValue',-70,'MarkerSize',3);
axis([0 maxDist/c -65 5]);
yticks([-60 -40 -20 0]);
ylabel('$h(t)$ [dB]');
xlabel('$t$');
xticks({});
ha = gca;
ha.YGrid = 'on';
hold on;

%% Animate propagation from center
dStep = 0.05; %m

for rad = 0:dStep:maxDist
    % radius is proportional to t
    
    subplot(4,1,1:3);
    delete(hr);
%     cen = PC.pos(1:2);
    for ii = 1:PC.n
        pos = [PC.pos(ii,1:2)-rad, rad*2, rad*2];
        hr(ii) = rectangle('Position',pos,'Curvature',[1,1],...
            'EdgeColor', [0 0 0 0.1+0.9*(1/(rad+1))], 'LineWidth', 2);
    end
    axis([-0.5*roomDims(1) 0.5*roomDims(1) -0.5*roomDims(2) 0.5*roomDims(2)]);
    
    subplot(4,1,4);
    delete(hs);
    sampEnd = ceil(rad*fs/c);
    hs = stem(tVec(1:sampEnd), db(rir(1:sampEnd)),'k','filled',...
        'BaseValue',-70,'MarkerSize',3);
    
%     pause(0.01);

    frame = getframe(gcf); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    if rad == 0 
        imwrite(imind,cm,gifname,'gif', 'Loopcount',inf,'DelayTime',tgif); 
    else 
        imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',tgif); 
    end
end

