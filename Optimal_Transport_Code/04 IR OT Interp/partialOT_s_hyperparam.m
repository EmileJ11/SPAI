`% Partial OT, with s as hyperparameter

% Author: Aaron Geldert
% Last modified: 18 May 2022

clear; clc; close all;
tic;
n = 50; % points in the distribution
dims = 3;
nVec = (1:n).';

% make dists with a virtual point
x = rand(n+1, dims); % assume dim 1 is mass, dims 2+ are coords
y = rand(n+1, dims);

% clear the virtual points (mass will be updated, pos is irrelevant)
x(n+1, :) = 0;
y(n+1, :) = 0;

% max amount of mass for partial transport
maxS = min(sum(x(1:n,1)), sum(y(1:n,1)));

xi = 0.01; % bounded scalar regularization parameter
q = 2; % distance exponent
C = zeros(n+1); % the extended cost matrix

sRatio = 0:0.05:1;
numS = numel(sRatio);
costs = zeros(numS, 1);

for iis = 1:numS
    s = sRatio(iis) * maxS;
    % virtual point mass
    x(n+1, 1) = sum(y(1:n,1)) - s; % untransported mass of other dist
    y(n+1, 1) = sum(x(1:n,1)) - s; 
    
    for ii = 1:n
        for jj = 1:n
            C(ii,jj) = norm(x(ii,2:dims) - y(jj,2:dims)).^q; % normal cost
        end
    end
    
    A = max(C,[],'all') + eps; % A > largest element of C
    C(n+1, :) = xi * ones(n+1, 1);
    C(:, n+1) = xi * ones(n+1, 1).'; % from eq. 1
    C(n+1, n+1) = 2*xi + A;
    
    cvx_begin quiet
    variable T(n+1,n+1) % coupling matrix
    expression d
        d = C .* T % extended cost * transport matrix
    minimize sum(d(:))
    subject to
        T >= 0; % matrix must be nonnegative
        sum(sum(T(1:n, 1:n))) == s; % total mass transported, w/o virtual pts
        T * ones(n+1,1) == x(:,1); % marginal sums kept, including virtual points
        T' *ones(n+1,1) == y(:,1);
    cvx_end
    
    costs(iis) = sum(sum(C.*T));

end
toc;

figure; plot(sRatio, costs);
xlabel('sRatio'); ylabel('C.*T');
title(['xi = ' num2str(xi)]);

%% Compare:

% With s = maximum, what is the amount of mass in Tin that is moved at a
% cost above xi?
% Does this mass correspond to the mass decrease from smax to sopt?

% Cex = C(1:end-1, 1:end-1);
% Tex = T(1:end-1, 1:end-1);
% Sex = sum(sum(Tex(Cex>xi)));
% Sopt = (maxS - Sex) / maxS

% smax = 24.94

% Yeah, this doesn't work... the assignment of masses to preserve marginals
% doesnt mean that all excess mass here is going to actually move to
% virtuals. It might move to other masses depending on s...


%% NOTES

% note that xi * (maxS - s) is the cost due to virtual points -
% deterministic

% note that the sum of T also scales depending on s...

% consider C.*T as 

% make a 2D surface of xi and s, normalized somehow

% You can see how much in C(orig) is at a higher cost than xi, because the
% budget for mapping s mass is not satisfied

% 
% Entropic regularization - can it remove artifacts or something?


