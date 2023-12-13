function PC = imageSourcePointCloud(roomDims, srcPos, maxOrder, doAtten)
% Returns the positions and magnitudes of image source in a cuboid room
% Assumes that origin [0 0 0] is in the center of room!
%   roomDims: [L W H] dimensions (m)
% Author: Aaron Geldert
% Last modified: 22 Aug 2022

if nargin < 4
    doAtten = 1;
end

% Calculate IS indices list
maxNum = 6 * 5^(maxOrder-1);
isInds = zeros(maxNum, 3);

% scan all possible unique indices and store
num = 0;
for x = -maxOrder:maxOrder
    k1 = maxOrder - abs(x);
    for y = -k1:k1
        k2 = k1 - abs(y);
        for z = -k2:k2
            num = num+1;
            isInds(num,:) = [x y z];
        end
    end
end

% remove blank rows
isInds = isInds(1:num,:);

% calculations
isOrders = vecnorm(isInds, 1, 2); % 1-norm is total order
PC.pos = srcPos.*(-1).^isInds + roomDims.*isInds;

PC.mass = ones(num,1).*(0.7071.^isOrders); % -3 dB per reflection
if doAtten ~= 0
    PC.mass = PC.mass ./ (vecnorm(PC.pos,2,2)+0.1); % 1/r attenuation with a regularization
end

PC.n = numel(PC.mass);

end

