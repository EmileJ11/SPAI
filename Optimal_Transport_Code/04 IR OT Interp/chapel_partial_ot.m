%% Partial Optimal Transport - Chapel method
% Reference: Chapel et al, Partial OT (2020)

% Author: Aaron Geldert
% Last modified: 18 May 2022

clear; clc; close all;

n = 4; % points in the distribution
dims = 3;
nVec = (1:n).';

% make dists with a virtual point
x = rand(n+1, dims); % assume dim 1 is mass, dims 2+ are coords
y = rand(n+1, dims);

% clear the virtual points (mass will be updated, pos is irrelevant)
x(n+1, :) = 0;
y(n+1, :) = 0;

% Don't normalize in this case!
% % normalize masses
% x(1:n,1) = x(1:n,1)./sum(x(1:n,1));
% y(1:n,1) = y(1:n,1)./sum(y(1:n,1));

% max amount of mass for partial transport
maxS = min(sum(x(1:n,1)), sum(y(1:n,1)));

% choose that 80% mass will be transported
sRatio = 0.8;
s = sRatio * maxS;

% virtual point mass
x(n+1, 1) = sum(y(1:n,1)) - s; % untransported mass of other dist
y(n+1, 1) = sum(x(1:n,1)) - s; 

%% Cost matrix
C = zeros(n+1); % the extended cost matrix

q = 2; % distance exponent
for ii = 1:n
    for jj = 1:n
        C(ii,jj) = norm(x(ii,2:dims) - y(jj,2:dims)).^q; % normal cost
    end
end

xi = mean(C,'all'); % a bounded scalar parameter - mean of C

A = max(C,[],'all') + 1; % A is greater than largest element of C
C(n+1, :) = xi * ones(n+1, 1);
C(:, n+1) = xi * ones(n+1, 1).'; % from eq. 1
C(n+1, n+1) = 2*xi + A;

% plots
xlims = [0.5 n+1.5];
figure;
subplot(223);
stem(x(1:n,1)); view([90 90]); grid on; hold on;
stem(n+1, x(end,1), 'md'); 
xlim(xlims); title(['Dist X: mass = ' num2str(sum(x(1:n,1)))]);

subplot(222);
stem(y(1:n,1)); grid on; hold on;
stem(n+1, y(end,1), 'md'); 
xlim(xlims); title(['Dist Y: mass = ' num2str(sum(y(1:n,1)))]);

subplot(224);
imagesc(C); title('Extended Cost Matrix');
xlabel('Y Indices'); ylabel('X Indices');
% set(gca, 'YDir', 'normal');

%% Optimization - s is input parameter
cvx_begin
    variable T(n+1,n+1) % coupling matrix
    expression d
        d = C .* T % extended cost * transport matrix
    minimize sum(d(:))
    subject to
        
        T >= 0; % matrix must be nonnegative
        
        sum(sum(T(1:n, 1:n))) == s; % total mass transported, w/o virtual pts
        
%         T(1:n, 1:n) * ones(n,1) <= x(1:n, 1); % constraints from initial PI definition
%         T(1:n, 1:n)' *ones(n,1) <= y(1:n, 1);
        
        T * ones(n+1,1) == x(:,1); % marginal sums kept, including virtual points
        T' *ones(n+1,1) == y(:,1);
        
cvx_end

%% results
% T = round(T);

subplot(221);
imagesc(T);
title('Transport Matrix');
xlabel('Y Indices'); ylabel('X Indices');


%% plot the 3D connections

T(T < 1e-6) = 0; % remove the unusably low values

figure;
scale = 200;
offset = 0.009;
scatter(x(1:n,2), x(1:n,3), (.001+x(1:n,1))*scale, 'b', 'filled');
hold on;
scatter(y(1:n,2), y(1:n,3), (.001+y(1:n,1))*scale, 'r', 'filled');

text(x(1:n, 2)+offset, x(1:n, 3), num2str((1:n)'), 'color', 'b');
text(y(1:n, 2)+offset, y(1:n, 3), num2str((1:n)'), 'color', 'r');

[row,col,val] = find(T);
for ii = 1:length(val)
    line([x(row(ii),2), y(col(ii),2)], [x(row(ii),3), y(col(ii),3)], 'Color','black');
end
xlabel('x'); ylabel('y');
title(['Transport ratio = ' num2str(sRatio)]);


