%% Arni SDM data 
% Created by: Aaron Geldert and Nils Meyer-Kahlen
% Last modified: 25 May 2022

filename = 'SDMData_arni_closed_r1.sofa';
D = SOFAload(filename);
fs = D.Data.SamplingRate; 
c = 343;

% 5 source positions
irs = squeeze(D.Data.IR).'; % samples x positions
doa = D.EmitterPosition; % samples x doa x positions

