function h = PCmapPlot(interpPCOT, PCgt, fs)
% Created by: Aaron Geldert & Nils Meyer-Kahlen
% Last modified: 25 Oct 2022
% plots several PCs on 3D map with an omni RIR plot
% 
% h = figure;

aziGrid = [0 30  60 90 120 150 180 -30 -60  -120 -150 -90 179.9]'/180*pi;
eleGrid = [-60 -30 0 30 60 ]'/180*pi;
linegray = 0.7;

subplot(6, 1, 1:3);
for k=1:length(aziGrid)
    ele=linspace(-pi/2, pi/2, 50);
    azi=ones(size(ele))*aziGrid(k);
    xyGrid = hammerAidhofProjection([azi', ele']);
    plot(xyGrid(:, 1), xyGrid(:, 2),'--', 'color', [1 1 1]*linegray);
    hold on
    text(xyGrid(25, 1),0, num2str(round(180/pi*aziGrid(k))), ...
        'fontsize', 8, 'color', [1 1 1]*linegray );
end

for k=1:length(eleGrid)
    azi=linspace(-pi, pi-.01, 50);
    ele=ones(size(azi))*eleGrid(k);
    xyGrid = hammerAidhofProjection([azi', ele']);
    plot(xyGrid(:, 1), xyGrid(:, 2),'--', 'color', [1 1 1]*linegray);
    hold on
    text(0, xyGrid(25, 2), num2str(180/pi*eleGrid(k)), ...
        'fontsize', 8, 'color', [1 1 1]*linegray );
end

PC1 = interpPCOT.PC1;
PC2 = interpPCOT.PC2;
PCk = interpPCOT.PCk;

scale = 1000; 
PC1sph = deg2rad(C2S(PC1.pos));
PC1hap = hammerAidhofProjection(PC1sph(:,1:2));
scatter(PC1hap(:,1), PC1hap(:,2), scale * PC1.mass, 'b.');

PC2sph = deg2rad(C2S(PC2.pos));
PC2hap = hammerAidhofProjection(PC2sph(:,1:2));
scatter(PC2hap(:,1), PC2hap(:,2), scale * PC2.mass, 'r.');

PCksph = deg2rad(C2S(PCk.pos));
PCkhap = hammerAidhofProjection(PCksph(:,1:2));
scatter(PCkhap(:,1), PCkhap(:,2), scale * PCk.mass, 'm.');

title('Image Sources');
xlabel('Azimuth'); xticks([]);
ylabel('Elevation'); yticks([]);

% IR plot down here
subplot(614);
maxT = 0.6;
[ir1, tvec] = monoIRfromPC(PC1, maxT, fs); ir1(ir1==0) = NaN;
[ir2, ~] = monoIRfromPC(PC2, maxT, fs); ir2(ir2==0) = NaN;
[irk, ~] = monoIRfromPC(PCk, maxT, fs); irk(irk==0) = NaN;
lastSample = max(find(ir1>0, 1, 'last'), find(ir2>0, 1, 'last'));
plotEndX = tvec(lastSample)*1.1e3;
stem(1e3*tvec, ir1, 'b','filled'); axis([0 plotEndX 0 .5]); ylabel('Pos 1');
subplot(615);
stem(1e3*tvec, irk, 'm','filled'); axis([0 plotEndX 0 .5]); ylabel('Interpolated');
subplot(616); 
stem(1e3*tvec, ir2, 'r','filled'); axis([0 plotEndX 0 .5]); ylabel('Pos 2');
xlabel('Time (msec)');

end

function [xy] = hammerAidhofProjection(aziEle)
    % Hammer-Aidhof projection
    azi = aziEle(:, 1);
    ele = aziEle(:, 2);
    hap = @(azi, ele) [-cos(ele).*sin(azi/2) 0.5 * sin(ele)] ./ ...
        (sqrt(1+cos(ele).*cos(azi/2)));
    xy = hap (mod (azi+ pi, 2 * pi) - pi, ele);
end