%% Partial Optimal Transport examples

% Method:
% 1. Consider all the 1-1 mappings permutations possible
% 2. Determine the cost of each:
%   - Distance of mapping
%   - Mass imbalance of mapping
%       - Multiply distance by ratio of masses?
% 3. Build up a partial permutation matrix until enough mass of each
% distribution has been assigned (minimizing cost)

% Author: Aaron Geldert
% Last modified: 10 May 2022

clear; clc; close all;

n = 5;
dims = 3;
nVec = (1:n).';

% make dists
x = rand(n, dims); % assume dim 1 is mass, dims 2+ are coords
y = rand(n, dims);

% normalize
x(:,1) = x(:,1)./sum(x(:,1));
y(:,1) = y(:,1)./sum(y(:,1));

%% cost matrix
costBal = 0.5; % 0:1, 0 = only distance, 1 = only mass
Cd = zeros(n); % distance cost
Cm = zeros(n); % mass imbalance cost
q = 2;
for ii = 1:n
    for jj = 1:n
        Cd(ii,jj) = norm(x(ii,2:dims) - y(jj,2:dims)).^q;
        Cm(ii,jj) = max(x(ii,1)./y(jj,1), y(jj,1)./x(ii,1)) - 1; % set so min is 0
    end
end
C = (Cd*(1-costBal) +1) .* (Cm*costBal +1);

% plots
figure;
subplot(231);
stem(x(:,1));

subplot(232);
stem(y(:,1));

subplot(234);
imagesc(Cd); title('Cd');
xlabel('Y Indices'); ylabel('X Indices');
set(gca, 'YDir', 'normal');

subplot(235);
imagesc(Cm); title('Cm');
xlabel('Y Indices'); ylabel('X Indices');
set(gca, 'YDir', 'normal');

subplot(236);
imagesc(C); title('C');
xlabel('Y Indices'); ylabel('X Indices');
set(gca, 'YDir', 'normal');

%% Optimization
m_partial = 0.7;
cvx_begin
    variables S(n,n) % coupling matrix
    expression d
    d = C .* S
    minimize sum(d(:))
    subject to

    % map points, disregard masses
    sum(S(:)) >= n*m_partial % get enough total mass assigned!
    S * ones(n,1) <= 1
    S' *ones(n,1) <= 1
    S >= 0

    % preserve mass
%     sum(S(:)) == 1
%     S * ones(n,1) <= x(:,1)
%     S' *ones(n,1) <= y(:,1)
%     S >= 0
    
cvx_end

S = round(S);

S(S < 1e-6) = 0;

subplot(233);
imagesc(S);
title('S');
xlabel('Y Indices'); ylabel('X Indices');
set(gca, 'YDir', 'normal');

%% plot the 3D connections
figure;
scale = 400;
offset = 0.009;
scatter(x(:,2), x(:,3), (.001+x(:,1))*scale, 'b', 'filled');
hold on;
scatter(y(:,2), y(:,3), (.001+y(:,1))*scale, 'r', 'filled');

text(x(:, 2)+offset, x(:, 3), num2str((1:n)'), 'color', 'b');
text(y(:, 2)+offset, y(:, 3), num2str((1:n)'), 'color', 'r');

[row,col,val] = find(S);
for ii = 1:length(val)
    line([x(row(ii),2), y(col(ii),2)], [x(row(ii),3), y(col(ii),3)], 'Color','black');
end
xlabel('x'); ylabel('y');
title(['m = ' num2str(m_partial)]);


