function xyz = randomPointOnSphere(ptNum, doPlot)
%RANDOMPOINTONTSPHERE Randomly selects a 3D unit vector
%   returns a cartesian coordinate [X, Y, Z], vector length is 1
% Created by: Aaron Geldert
% Last modified: 15 Oct 2022

    if nargin < 1
        ptNum = 1;
    end
    if nargin < 2
        doPlot = 0;
    end
    
    xyz = randn(ptNum, 3);
    xyz = xyz./vecnorm(xyz,2,2);
    
    if doPlot ~= 0
        scatter3(xyz(:,1), xyz(:,2), xyz(:,3), '.');
        axis equal;
    end
end