% Pulse signal example for the thesis report.

% Add path
paths = genpath('C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Matlab\');
addpath(paths);

% Options for generation of images
imagGenFlag = true; % Generate images

% Parameters
pulseCoefMat(:, :, 1) = ...
    [1 -1;
    0 0;
    0 0;
    -1 1];

pulseCoefMat(:, :, 2) = ...
    [0 0;
    1 0;
    0 1i;
    0 0];

pulseLimits = 20*...
    [1 2;
    3 4;
    5 7;
    9 10];
freqs = [2; 4];
channels = [2 1];
sampleRate = 10000;

% Generate signal
[signal, outOfRange] = pulseCoefMat2signal(pulseCoefMat, pulseLimits, freqs, sampleRate, 0, 250, 'type_marker', 'time', 'type_pulseLimits', 'time');
numSamples = size(signal, 1);
t = (0:numSamples-1)/sampleRate;

% Represent signal in time
fig = figure;
ax = axes(fig);
plot(ax, t, signal)
ax.XTick = sort(pulseLimits(:));
ax.XLabel.String = 'Time (s)';

numChannels = numel(channels);
lengendStrings = cell(numChannels, 1);
for c = 1:numChannels
    lengendStrings{c} = sprintf('Channel %d', channels(c));
end
legend(lengendStrings);

if imagGenFlag
    name = 'pulseSignal';
    print(fig, ['C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\', name, '.eps'], '-depsc');
    print(fig, ['C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\', name, '.emf'], '-dmeta');
end

% Perform DFT
X = fft(signal(:,2));
duration = numSamples/sampleRate;
df = 1/duration;
f = (0:numSamples - 1)*df;

f_max = sampleRate;
f_2 = f_max/2;
f_shift = fftshift(mod(f + f_2, f_max) - f_2);
X_shift = fftshift(X);

fig = figure;
ax = axes(fig);
plot(ax, f_shift, abs(X_shift));
ax.XLim = [-max(freqs)*1.3, max(freqs)*1.3];
ax.YTick = [];
ax.XTick = sort([-freqs; 0; freqs]);
ax.XLabel.String = 'Frequency (Hz)';

if imagGenFlag
    name = 'DFT_channel2';
    print(fig, ['C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\', name, '.eps'], '-depsc');
    print(fig, ['C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\', name, '.pdf'], '-dpdf');
    print(fig, ['C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\', name, '.emf'], '-dmeta');
end

% Apply filter. Draw Mask
f_filter = 2;
width = 1;

mask = zeros(numSamples, 1);
mask(f_shift > f_filter - width/2 & f_shift < f_filter + width/2) = max(abs(X_shift))*1.3;

ax.NextPlot = 'Add';
plot(ax, f_shift, mask)

if imagGenFlag
    name = 'DFT_channel2_filter';
    print(fig, ['C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\', name, '.eps'], '-depsc');
    print(fig, ['C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\', name, '.pdf'], '-dpdf');
    print(fig, ['C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\', name, '.emf'], '-dmeta');
end

% IQ Signal
% Filter
pitchFreqSamp = f_filter*numSamples/sampleRate;
filterWidthSamp = width*numSamples/sampleRate;
lowFreq = floor(pitchFreqSamp - filterWidthSamp);
highFreq = ceil(pitchFreqSamp + filterWidthSamp);
X(1:lowFreq-1, :) = 0;
X(highFreq+1:end, :) = 0;
x_filtered = ifft(X)*2;

% Downconvert by the specified frequency
A = exp(-1i*2*pi*f_filter/sampleRate*(0:numSamples-1)');
iq = x_filtered.*A;

fig = figure;
ax = axes(fig);
plot(ax, t, real(iq), t, imag(iq));
ax.XLabel.String = 'Time (s)';

if imagGenFlag
    name = 'IQ';
    print(fig, ['C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\', name, '.eps'], '-depsc');
    print(fig, ['C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\', name, '.pdf'], '-dpdf');
    print(fig, ['C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\', name, '.emf'], '-dmeta');
end

% Find the temporal shift with a correlation

% Generate mask from limits
numPulses = size(pulseLimits, 1);
mask = zeros(numSamples, 1);
for p = 1:numPulses
    mask(t >= pulseLimits(p, 1) & t <= pulseLimits(p, 2)) = 1;
end
shiftTime = 2;
shiftSamp = floor(shiftTime*sampleRate);
mask = shiftAndCrop(mask, shiftSamp, numSamples);

fig = figure;
ax = axes(fig);
plot(ax, t, abs(iq), t, mask)

if imagGenFlag
    name = 'correlation';
    print(fig, ['C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\', name, '.eps'], '-depsc');
    print(fig, ['C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\', name, '.pdf'], '-dpdf');
    print(fig, ['C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\', name, '.emf'], '-dmeta');
end



