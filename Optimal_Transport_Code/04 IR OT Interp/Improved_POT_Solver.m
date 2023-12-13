%% A partial OT solution that does not require iterative solutions!
% Test uses 2D point distributions
% Created by: Aaron Geldert
% Last modified: 7 Dec 2022

% Wow, speed is improved! Fewer parameters too!

clear; clc; close all;
n = 5;

PC1.dist = 5+randn(n,1);
PC1.dir = 360*rand(n,1);
PC1.mass = 0.1+0.9*rand(n,1);

PC2.dist = 2.5*rand(n,1);
PC2.dir = 360*randn(n,1);
PC2.mass = 0.1+0.9*rand(n,1);

PC1.pos = S2C([PC1.dir, zeros(n,1), PC1.dist]);
PC2.pos = S2C([PC2.dir, zeros(n,1), PC2.dist]);

% cvx definitions
C = zeros(n);
for ii = 1:n
    for jj = 1:n
        C(ii,jj) = norm(PC1.pos(ii,:) - PC2.pos(jj,:)).^2;
    end
end

dummyCost = 24;
Cx = [C dummyCost*ones(n,1)];
Cx = [Cx; dummyCost*ones(1,n+1)];
Cx(n+1,n+1) = max(max(Cx))+1;

tic;
cvx_begin
    cvx_solver sedumi
    variable Tx(n+1,n+1)
    expression d
        d = Cx .* Tx
    minimize sum(d(:))
    subject to
        sum(Tx(1:n,:),2) == PC1.mass
        sum(Tx(:,1:n),1).' == PC2.mass
        Tx >= 0
        Tx(n+1,n+1) == 0 % redundant
cvx_end
toc;

% disp(Cx);
Tx = full(Tx);
Txd = [[PC2.mass.', NaN]; Tx]; % add top row
Txd = [[NaN; PC1.mass; NaN], Txd]; % add left col

figure(1);
subplot(121); imagesc(Cx);
subplot(122); imagesc(Tx);

%% brute force permutation vector
map = zeros(n,1);
for ii = 1:n
    map(ii) = find(T(ii,:) > 0.5);
end
disp(map)

%% plot
gifname = 'ot2d_example2.gif';
num = 120;
tgif = 4/num;

figure(2); clf;
scatter(PC1.pos(:,1), PC1.pos(:,2),'MarkerEdgeColor',[.1 .1 .9]);
hold on;
scatter(PC2.pos(:,1), PC2.pos(:,2),'MarkerEdgeColor',[.9 .1 .1]);
xticks({}); yticks({});

% for ii = 1:n
%     line([PC1.pos(ii,1), PC2.pos(map(ii),1)],...
%         [PC1.pos(ii,2), PC2.pos(map(ii),2)],...
%         'Color', [.5 .5 .5 ]);
% end
hk = scatter(0,0);
for k = 0:0.01:1
    delete(hk);
    hk = scatter((1-k)*PC1.pos(:,1)+k*PC2.pos(:,1),...
        (1-k)*PC1.pos(:,2)+k*PC2.pos(:,2),'filled',...
        'MarkerFaceColor',[0.2+0.7*k 0.2 0.9-0.7*k]);
%     pause(0.02);
    frame = getframe(gcf); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    if k == 0 
        imwrite(imind,cm,gifname,'gif', 'Loopcount',inf,'DelayTime',tgif); 
    else 
        imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',tgif); 
    end
end
for k = 1:-0.01:0
    delete(hk);
    hk = scatter((1-k)*PC1.pos(:,1)+k*PC2.pos(:,1),...
        (1-k)*PC1.pos(:,2)+k*PC2.pos(:,2),'filled',...
        'MarkerFaceColor',[0.2+0.7*k 0.2 0.9-0.7*k]);
%     pause(0.02);
frame = getframe(gcf); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',tgif); 
end