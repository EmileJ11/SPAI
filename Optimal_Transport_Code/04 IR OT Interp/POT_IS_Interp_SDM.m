%% Partial Optimal Transport on SDM Data
% Author: Aaron Geldert
% Last modified: 13 May 2022

clear; clc; close all;

% Lets examine interpolation between positions 18:20
ir1 = load('ir_18.mat');
ir2 = load('ir_19.mat');
ir3 = load('ir_20.mat');

db2en = @(dbval) 10.^(dbval./10); % /10 for energy
isMag1 = db2en(ir1.masses);
isMag2 = db2en(ir3.masses);
isPos1 = ir1.positions;
isPos2 = ir3.positions;

% % do some alterations (add noise, etc.)
% magError = 0.3;
% isMag1 = isMag1 .* ((1-magError/2)+ magError*rand(size(isMag1)));
% isMag2 = isMag2 .* ((1-magError/2)+ magError*rand(size(isMag2)));
% 
% posError = 0.2;
% isPos1 = isPos1 .* ((1-posError/2)+ posError*rand(size(isPos1)));
% isPos2 = isPos2 .* ((1-posError/2)+ posError*rand(size(isPos2)));

% select n most prominent image sources as dists
n = 128;
dims = 4;
[~,xi] = sort(isMag1, 'descend');
xi = xi(1:n);
x = [isMag1(xi), isPos1(xi,:)]; % dim 1 is mass, dims 2:4 are coords

[~,yi] = sort(isMag2, 'descend');
yi = yi(1:n);
y = [isMag2(yi), isPos2(yi,:)];

% normalize
x(:,1) = x(:,1)./sum(x(:,1));
y(:,1) = y(:,1)./sum(y(:,1));

% cost matrix
costBal = 0.4; % 0:1, 0 = only distance, 1 = only mass
Cd = zeros(n); % distance cost
Cm = zeros(n); % mass imbalance cost
q = 2;
for ii = 1:n
    for jj = 1:n
        Cd(ii,jj) = norm(x(ii,2:dims) - y(jj,2:dims)).^q; % distance cost
        Cm(ii,jj) = max(x(ii,1)./y(jj,1), y(jj,1)./x(ii,1)) - 1; % set so min is 0
    end
end
C = (Cd*(1-costBal) +1) .* (Cm*costBal +1);

%% optimizations! 
% 1) Partial OT
m_partial = 0.7;
cvx_begin
    variables S1(n,n) % coupling matrix
    expression d
        d = C .* S1;
    minimize sum(d(:))
    subject to
        sum(S1(:)) >= n*m_partial; % get enough total mass assigned!
        S1 * ones(n,1) <= 1;
        S1' *ones(n,1) <= 1;
        S1 >= 0;
cvx_end
S1orig = S1; % original, rounding errors present
S1 = round(S1); % binary matrix 

%% plot
figure;
scale = 500;
offset = 0.005;
scatter3(x(:,2), x(:,3), x(:,4), x(:,1)*scale, 'b', 'filled');
hold on;
scatter3(y(:,2), y(:,3), y(:,4), y(:,1)*scale, 'r', 'filled');

[row,col,val] = find(S1);
for ii = 1:length(val)
    line([x(row(ii),2), y(col(ii),2)], [x(row(ii),3), y(col(ii),3)], [x(row(ii),4), y(col(ii),4)], 'Color','black');
end

xlabel('x'); ylabel('y'); zlabel('z');
title(['m = ' num2str(m_partial)]);





