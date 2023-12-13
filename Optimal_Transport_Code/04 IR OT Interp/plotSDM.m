function [outputArg1,outputArg2] = plotSDM(h, SDM, plotTitle, axlim, scale)
%plotSDM Summary of this function goes here
%   Detailed explanation goes here
% Created by: Aaron Geldert
% Last modified: 23 May 2022

scatter3(h, 0,0,0,'*'); % origin  (listener position)
hold on; grid on;
scatter3(h, SDM.pos(:,1), SDM.pos(:,2), SDM.pos(:,3),...
    1+SDM.mass * scale, 'ko','filled');
view(h, [0 90]);
xlabel(h, 'X'); ylabel(h, 'Y'); zlabel(h, 'Z');
title(h, plotTitle);
axis(axlim);

end

