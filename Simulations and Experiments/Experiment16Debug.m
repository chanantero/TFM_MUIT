%% 
dataPathName = 'Data/';
IDload = 'Lab_29-08-2018';
signalStructLoaded = load([dataPathName, 'signalStruct_', IDload, '.mat']);
freqStep = signalStructLoaded.freqs(2) - signalStructLoaded.freqs(1);
T = 1/freqStep;
numChannels = 96;
numSimultChan = signalStructLoaded.numSimultChan;
pulseLimits = signalStructLoaded.pulseLimits;
freqs = signalStructLoaded.freqs;
numFreqs = length(freqs);
coefMat = signalStructLoaded.coefMat;
filename = './Data\acousticPathReprodSignal_Lab_29-08-2018.wav';
info = audioinfo(filename);
sampleRate = info.SampleRate;

samplesPerFrame = sampleRate*2;
numFrames = ceil(info.TotalSamples/samplesPerFrame);

signProv = signalProvider;
signProv.FileName = filename;
signProv.SamplesPerFrame = samplesPerFrame; % The maximum an audioDeviceWriter allows

testSignalNS = zeros(numFrames*samplesPerFrame, 1);
ind = 1-samplesPerFrame:0;
for fr = 1:numFrames
    x_curr = step(signProv);
    ind = ind + samplesPerFrame;
    testSignalNS(ind) = x_curr(:, 1);
    fprintf('%d/%d\n', fr, numFrames)
end
release(signProv)

load('./Data/acousticPathRecordSignal_Lab_29-08-2018_corregido.mat');

recSignalTest = recSignal';

numSamp = length(testSignalNS);
t = (0:numSamp - 1)/sampleRate;
ax = axes(figure, 'NextPlot', 'Add');
plot(ax, t(t<12), testSignalNS(t<12))
plot(ax, t(t<12), recSignalTest((t<12), 1)')

freqsPlus = [freqs, freqs(end)+1:2000];

ind = t > pulseLimits(1,1) + 1 & t < pulseLimits(1,2) - 1;
fragmentX = testSignalNS(ind);
specX = freqz(fragmentX, 1, freqsPlus, sampleRate);
fragmentY = recSignalTest(ind, 1);
specY = freqz(fragmentY, 1, freqsPlus, sampleRate);

ax = axes(figure, 'NextPlot', 'Add');
plot(ax, freqsPlus, abs(specX))
plot(ax, freqsPlus, abs(specY))

coefMatFlag = permute(coefMat ~= 0, [3 2 1]);
flag = coefMatFlag(:, 1, 1);

Hofdm = zeros(length(freqs), 1);
Hofdm(flag) = specY(flag)./specX(flag);
Hofdm_script = squeeze(FR(1, 1, :));

ax = axes(figure, 'NextPlot', 'Add');
plot(ax, freqs, abs(Hofdm))
plot(ax, freqs, abs(Hofdm_script))

durSign = 4; % Duration of tone for time processing
t = (0:ceil(durSign*sampleRate)-1)/sampleRate;
NSsignal = chirp(t, min(freqs), durSign, max(freqs) + 50);
prefixDuration = 2;
prefixNumSamples = ceil(prefixDuration*sampleRate);
NSsignal = [zeros(1, prefixNumSamples), NSsignal, zeros(1, prefixNumSamples )];
t = (0:length(NSsignal)-1)/sampleRate;

load([dataPathName, 'recSignals_Lab_29-08-2018_corregido.mat'])

X = freqz(NSsignal, 1, freqsPlus, sampleRate);
Y = freqz(recSignalNS(1,:), 1, freqsPlus, sampleRate);
Hchirp = Y./X;

ax = axes(figure, 'NextPlot', 'Add');
plot(ax, freqsPlus, abs(X))
plot(ax, freqsPlus, abs(Y))

ax = axes(figure, 'NextPlot', 'Add');
plot(ax, freqs, abs(Hofdm))
plot(ax, freqs, abs(Hchirp))

ax = axes(figure, 'NextPlot', 'Add');
plot(ax, freqs(flag), rad2deg(unwrap(angle(Hofdm(flag)))))
plot(ax, freqs, rad2deg(unwrap(angle(Hchirp))))

phase = unwrap(angle(Hchirp));
delay = -phase./(2*pi*freqs);
ax = axes(figure);
plot(ax, freqs, delay)

freqDif = diff(freqs);
phaseDif = diff(phase);
deriv = phaseDif./freqDif;
groupDelay = -deriv/(2*pi);
ax = axes(figure);
plot(ax, freqs(1:end-1), groupDelay)

phaseOFDM = unwrap(angle(Hofdm(flag)));
delayOFDM = -phaseOFDM'./(2*pi*freqs(flag));
ax = axes(figure);
plot(ax, freqs(flag), delayOFDM)

freqDif = diff(freqs(flag));
phaseDif = diff(phaseOFDM);
deriv = phaseDif./freqDif;
groupDelayOFDM = -deriv/(2*pi);
ax = axes(figure);
freqsAux = freqs(flag);
plot(ax, freqsAux(1:end-1), groupDelayOFDM)

ax = axes(figure, 'NextPlot', 'Add');
plot(ax, NSsignal)
plot(ax, recSignalNS(1,:))

ax = axes(figure, 'NextPlot', 'Add');
plot(ax, fragmentX)
plot(ax, fragmentY)

%% Estimación OFDM paso a paso
ind = t > pulseLimits(1,1) + 1 & t < pulseLimits(1,2) - 1;
fragmentX = testSignalNS(ind);
specX = freqz(fragmentX, 1, freqs, sampleRate);
fragmentY = recSignalTest(ind, 1);
specY = freqz(fragmentY, 1, freqs, sampleRate);

coefMatFlag = permute(coefMat ~= 0, [3 2 1]);
flag = coefMatFlag(:, 1, 1);

Hofdm = zeros(length(freqs), 1);
Hofdm(flag) = specY(flag)./specX(flag);

%% Estimación