%% THE PROBLEM WITH UNBALANCED TRANSPORT - SPONTANEOUS MASS

pc1.pos = [1 0 0;...
           2 0 0];
pc1.mass = [2; 2];
pc1.n = 2;


pc2.pos = [5 0 0;...
           6 0 0];
pc2.mass = [1; 1];
pc2.n = 2;

xi = 5;
pcot = PCOT(pc1, pc2, xi^2);

[pctm, pcsm] = pcot.interpPC(0.5);

%%

pc1.pos = [1 0 0;...
           2 0 0;...
           5 0 0;...
           9 0 0];
pc1.mass = [6; 7; 4; 2];
pc1.n = 4;


pc2.pos = [2 0 0;...
           3 0 0;...
           5 0 0];
pc2.mass = [3; 8; 9];
pc2.n = 3;

xi = 2;
pcot = PCOT(pc1, pc2, xi^2);
[pctm, pcsm] = pcot.interpPC(0.5);

pcot.T = round(pcot.T);
pcot.Tx = round(pcot.Tx);
disp(pcot.Tx)

%% Some plots
close all;

figure(1);
subplot(211); stem(pc1.pos(:,1),pc1.mass, 'filled'); 
axis([0 10 -1 10]); title('Distribution 1'); grid on; ylabel('Mass');

subplot(212); stem(pc2.pos(:,1),pc2.mass , 'filled'); 
axis([0 10 -1 10]); title('Distribution 2'); grid on; xlabel('Position'); ylabel('Mass');

figure(2);
imagesc(pcot.Cx); 
title('Extended Cost Matrix');
xlabel('Dist 2 Indices'); ylabel('Dist 1 Indices'); 
xticks(1:4); yticks(1:5); 

ytl = cell(pc1.n+1, 1);
for yi = 1:pc1.n
    ytl{yi} = [num2str(yi) ': (' num2str(pc1.pos(yi)) ',', num2str(pc1.mass(yi)) ')'];
end
ytl{yi+1} = 'Dummy';

xtl = cell(pc2.n+1, 1);
for xi = 1:pc2.n
    xtl{xi} = [num2str(xi) ': (' num2str(pc2.pos(xi)) ',', num2str(pc2.mass(xi)) ')'];
end
xtl{xi+1} = 'Dummy';

xticklabels(xtl); yticklabels(ytl);
set(gca,'xaxisLocation','top');
colorbar;
set(gca,'ColorScale','log');

for ii=1:pc2.n+1
    for jj=1:pc1.n+1
        text(ii, jj, num2str(pcot.Cx(jj,ii)));
    end
end

figure(3);
imagesc(pcot.Tx); 
title('Extended Transport Matrix');
xlabel('Dist 2 Indices'); ylabel('Dist 1 Indices'); 
xticks(1:4); yticks(1:5); 

ytl = cell(pc1.n+1, 1);
for yi = 1:pc1.n
    ytl{yi} = [num2str(yi) ': (' num2str(pc1.pos(yi)) ',', num2str(pc1.mass(yi)) ')'];
end
ytl{yi+1} = 'Dummy';

xtl = cell(pc2.n+1, 1);
for xi = 1:pc2.n
    xtl{xi} = [num2str(xi) ': (' num2str(pc2.pos(xi)) ',', num2str(pc2.mass(xi)) ')'];
end
xtl{xi+1} = 'Dummy';

xticklabels(xtl); yticklabels(ytl);
set(gca,'xaxisLocation','top');

for ii=1:pc2.n+1
    for jj=1:pc1.n+1
        if pcot.Tx(jj,ii) >0
            text(ii, jj, num2str(pcot.Tx(jj,ii)));
        end
    end
end
