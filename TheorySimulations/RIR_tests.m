%% Study the difference between the ideal IR and the one returned by RIR-generator

%%

numSamp = 4096;
fs = 44100;
distance = 5;

c = 340;
atten = 1/distance;
delay = distance/c;
t = (0:numSamp-1)/fs;

IR_ideal = zeros(1, numSamp);
ind = floor(delay*fs);
IR_ideal(ind) = atten;

% Generator
offset = [5 5 5];
r = [distance 0 0] + offset;
s = [0 0 0] + offset;
roomDim = [10 10 10];
IR_rir = rir_generator(c, fs, r, s, roomDim, 0, numSamp);

% % Plot both IR
% ax = axes(figure);
% plot(ax, t, IR_ideal, t, IR_rir)

% Peform spectral analysis
f = (0:numSamp - 1)*(fs/numSamp);
FR_ideal = fft(IR_ideal);
FR_rir = fft(IR_rir);

% axFR = axes(figure);
% plot(axFR, f, abs(FR_ideal), f, abs(FR_rir))
% yyaxis right
% plot(axFR, f, rad2deg(angle(FR_ideal)), '--b', f, rad2deg(angle(FR_rir)), '--r')
% axFR.YAxis(1).Limits = [0 1.5];

axPhase = axes(figure);
plot(axPhase, f, rad2deg(unwrap(angle(FR_ideal))), f, rad2deg(unwrap(angle(FR_rir))))

rel = FR_rir./FR_ideal;
relPhase = angle(rel);
plot(axPhase, f, rad2deg(unwrap(relPhase)))

% df = fs/numSamp;
% groupDelay_ideal = diff(unwrap(angle(FR_ideal)))/(2*pi*df);
% groupDelay_rir = diff(unwrap(angle(FR_rir)))/(2*pi*df);
% 
% axGroupDelay = axes(figure);
% plot(axGroupDelay, f(1:end-1), groupDelay_ideal, f(1:end-1), groupDelay_rir)



%% Different distances

distance = (1:1:100)';
numDist = numel(distance);
numSamp = 20000;
fs = 44100;

c = 340;
atten = 1./distance;
delay = distance/c;
t = (0:numSamp-1)/fs;

IR_ideal = zeros(numDist, numSamp);
ind = floor(delay*fs);
for d = 1:numDist
    IR_ideal(d, ind(d)) = atten(d);
end

% Generator
offset = [5 5 5];
r = [distance zeros(numDist, 2)] + repmat(offset, numDist, 1);
s = [0 0 0] + offset;
roomDim = [10 10 10];
IR_rir = rir_generator(c, fs, r, s, roomDim, 0, numSamp);

df = fs/numSamp;
f = (0:numSamp - 1)*df;
FR_ideal = fft(IR_ideal');
FR_rir = fft(IR_rir');

groupDelay_ideal = diff(unwrap(angle(FR_ideal)), 1)/(2*pi*df);
groupDelay_rir = diff(unwrap(angle(FR_rir)), 1)/(2*pi*df);

selFreqs = [440, 600, 800];
ind = zeros(numel(selFreqs), 1);
for k = 1:numel(selFreqs)
    ind(k) = find(f <= selFreqs(k), 1, 'last'); % closest indices to set of frequencies
end

% ax = axes(figure);
% plot(ax, distance, groupDelay_ideal(ind, :)', '--b', distance, groupDelay_rir(ind, :)')

rel = FR_rir./FR_ideal;
relPhase = angle(rel);
relAbs = abs(rel);

% ax2 = axes(figure);
% plot(ax2, distance, rad2deg(relPhase(ind, :))')

% ax3 = axes(figure);
% plot(ax3, distance, rad2deg(unwrap(angle(FR_ideal(ind, :))', [], 1)), distance, rad2deg(unwrap(angle(FR_rir(ind, :))', [], 1)))

ax4 = axes(figure);
plot(ax4, f, rad2deg(unwrap(relPhase)))
ax4.XLim = [0 1000];

%% Understanding the reverberation time and the reflexion coefficient

h = 1; d = 2;
sourcePos = [1 1 h];
recPos = [sourePos(1) + d 1 h];
c = 340;
fs = 44100;
roomDim = [300 200 2*h];
Beta = [0 0 0 0 1 1];
numSampIR = 4*1024;
t = (0:numSampIR-1)/fs;

h = rir_generator(c, fs, recPos, sourcePos, roomDim, Beta, numSampIR);
ax = axes(figure);
plot(ax, t, h)

hmod = h; hmod(hmod<0) = 0;
[maxtab, mintab] = peakdet(hmod, 0.01, t);

numRebot = 1:4;
distReb = sqrt((2*h*numRebot).^2 + d^2);
distReb/c


[h, betaHat] = rir_generator(c, fs, recPos, sourcePos, roomDim, 0.16, numSampIR);
h2 = rir_generator(c, fs, recPos, sourcePos, roomDim, betaHat*ones(1,6), numSampIR);
ax = axes(figure);
plot(ax, t, h, t, h2)
