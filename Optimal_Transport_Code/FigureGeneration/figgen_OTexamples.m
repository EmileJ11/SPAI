%% Nice OT examples
clc;
x = 1:1000;

p1 = 2*[hann(500); zeros(500,1)];
p1 = p1 - [zeros(50,1); hann(350); zeros(600,1)];

p2 = [zeros(500,1); hann(500);];
p2 = p2 + [zeros(550,1); hann(250).^2; zeros(200,1)];
p2 = p2 - 0.5*[zeros(580,1); flattopwin(400).^2; zeros(20,1)];

p1 = 0.5*p1;
p2 = 0.5*p2;

area(x,p1,'FaceColor',[0.7 0.05 0.1],'FaceAlpha',0.2);
hold on;
area(x,p2,'FaceColor',[0.1 0.05 0.7],'FaceAlpha',0.2);

OT = OptimalTransport1D(p1,p2);
[B,A] = butter(8,0.1);

gifname = 'ot1d_example1.gif';
num = 120;
tgif = 4/num;

figure(1);
xticks({}); yticks({});
legend('Distribution 1','Distribution 2','Location','northwest');
set(gcf,'Position',[100,100,400,200]);

for k = 0:0.01:1
    delete(hk);
    pk = OT.interpolate(k);
    pk = filtfilt(B,A,pk);
    pk(pk<0) = 0;
    hk = area(x,pk,'FaceColor',[0.7-0.6*k 0.05 0.1+0.6*k],'FaceAlpha',0.5,'DisplayName','Transported');
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
    pk = OT.interpolate(k);
    pk = filtfilt(B,A,pk);
    pk(pk<0) = 0;
    hk = area(x,pk,'FaceColor',[0.7-0.6*k 0.05 0.1+0.6*k],'FaceAlpha',0.5,'DisplayName','Transported');
%     pause(0.02);
    frame = getframe(gcf); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',tgif); 
end

%% Point clouds

clear; clc; close all;
n = 64;

PC1.dist = 5+randn(n,1);
PC1.dir = 360*rand(n,1);

PC2.dist = 2.5*rand(n,1);
PC2.dir = 360*randn(n,1);

PC1.pos = S2C([PC1.dir, zeros(n,1), PC1.dist]);
PC2.pos = S2C([PC2.dir, zeros(n,1), PC2.dist]);

% cvx beauty

C = zeros(n);
for ii = 1:n
    for jj = 1:n
        C(ii,jj) = norm(PC1.pos(ii,:) - PC2.pos(jj,:)).^2;
    end
end

cvx_begin quiet
    cvx_solver sedumi
    variable T(n,n)
    expression d
        d = C .* T
    minimize sum(d(:))
    subject to
        T * ones(n,1) == ones(n,1)
        T' *ones(n,1) == ones(n,1)
        T >= 0
        T <= 1
cvx_end

% brute force permutation vector
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


