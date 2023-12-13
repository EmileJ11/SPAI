% Lets look at these results:
% Author: Aaron Geldert
% Last modified: 5 Oct 2022

close all;
clear;

pathout = '/Users/aarongeldert/Documents/MATLAB/THESIS/Exported Results/IS_movingRcv_staticSrc/';

figure(1);
for pairSep = 2:20
    figure(pairSep-1);
    sgtitle(['Pair Separation: ' num2str(pairSep)]);
    filename = [pathout 'IS_movingRcv_staticSrc_sep' num2str(pairSep) '.mat'];
    LD(pairSep-1) = load(filename);
    subplot(121);
    boxplot(LD(pairSep-1).linError.');
    ylim([0 1]);
    title('Linear Interp', 'Units', 'normalized', 'Position', [0.5, 0.92, 0]);
    subplot(122);
    boxplot(LD(pairSep-1).potError.');
    ylim([0 1]);
    title('OT Interp', 'Units', 'normalized', 'Position', [0.5, 0.92, 0]);
end