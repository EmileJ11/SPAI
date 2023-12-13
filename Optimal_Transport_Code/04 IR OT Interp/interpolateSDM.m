function SDMint = interpolateSDM(SDM1, SDM2, k, xi, sRatio)
%INTERPOLATESDM Interpolates between two SDM point clouds at k
% Uses Chapel's partial OT method
% Created by: Aaron Geldert
% Last modified: 24 May 2022

if nargin < 3
    k = 0.5;
end
if nargin < 4
    xi = 0.05;
end
if nargin < 5
    sRatio = 0.75;
end

% Add virtual points
m1 = [SDM1.mass; 0];
m2 = [SDM2.mass; 0];

% Determine max amount of mass (s)
maxS = min(sum(m1), sum(m2));
s = sRatio * maxS;

m1(end) = sum(SDM2.mass) - s;
m2(end) = sum(SDM1.mass) - s;

% Prepare cost matrix (distance-based)
q = 2;
Corig = zeros(SDM1.n, SDM2.n);
for ii = 1:SDM1.n
    for jj = 1:SDM2.n
        Corig(ii,jj) = norm(SDM1.pos(ii,:) - SDM2.pos(jj,:)).^q; % normal cost
    end
end

A = max(Corig,[],'all') + 1; % A > largest element of C

xi = mean(Corig,'all'); % overwrite xi parameter

Cv1 = xi * ones(SDM1.n, 1); 
Ctemp = horzcat(Corig, Cv1); % add virtual point costs

Cv2 = xi * ones(SDM2.n, 1).'; 
C = vertcat(Ctemp, horzcat(Cv2, A));

% Optimization
cvx_begin
    variable T(SDM1.n+1, SDM2.n+1)
    expression d
        d = C .* T
    minimize sum(d(:))
    subject to
        T >= 0; % matrix must be nonnegative
        sum(sum(T(1:end-1, 1:end-1))) == s; % total mass transported, w/o virtual pts
        T * ones(SDM2.n+1,1) == m1; % marginal sums kept, including virtual points
        T' *ones(SDM1.n+1,1) == m2;
cvx_end

T(T < 1e-5) = 0; % round down to 0

% Interpolate using transport plan
[DOAint, Mint] = interpSDM(k, SDM1, SDM2, T, 1);

SDMint.fs = SDM1.fs;
SDMint.tVec = SDM1.tVec;
SDMint.pos = DOAint;
SDMint.dist = sqrt(DOAint(:,1).^2 + DOAint(:,2).^2 + DOAint(:,3).^2);
SDMint.mass = Mint;
SDMint.n = numel(Mint);

end

