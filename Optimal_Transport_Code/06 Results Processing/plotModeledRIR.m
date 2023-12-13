function h = plotModeledRIR(PC, fs, stemColor)
% DEPRECATED

assert(nargin <= 2, 'Not enough input arguments');
assert(isstruct(PC),'PC must be a struct');
assert(isscalar(fs),'fs must be a scalar');
% assert(isvalid(h),'h must be a handle to Axes');

sf = hann(round(0.004 * fs));
[rir, tvec] = monoIRfromPC(PC, 0.2, fs);
edc = conv(rir, sf, 'same');

plotTime = 0.04; % ms of plot desired
dbGap = 40;
dbBln = 35; linBln = 10.^(dbBln/20);
dbOffset = 3;
dbTickGap = -20;
axisLim1 = [0 1e3*plotTime -dbGap-dbBln 0];
axisLim2 = [0 1e3*plotTime -3*dbGap-dbBln 0];
stemWidth = 1.5;

dbRir = db(rir) + dbOffset;
dbRir(dbRir < -dbBln) = -Inf;
dbEdc = db(edc) + dbOffset;
dbEdc(dbEdc < -dbBln) = -Inf;

greyColor = [.75 .75 .75];

stem(1e3*tvec, dbRir,'filled','k',...
    'BaseValue', -dbBln, 'MarkerSize', 3, 'LineWidth', stemWidth);
hold on;
plot(1e3*tvec, dbEdc,'k');
axis(axisLim1); yticklabels({}); xticklabels({});
h = gca;
set(h,'YTick',[])


end