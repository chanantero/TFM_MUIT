%% Windows comparison

% Parameters
pulseDuration = 1;
fftDuration = 20;
sampleRate = 1024;

% Derived parameters
numPulseSamples = round(pulseDuration * sampleRate);
fftSize = round(fftDuration*pulseDuration*sampleRate);

deltaT = 1/sampleRate;
t = (0:numPulseSamples-1)/sampleRate;

deltaF = sampleRate/fftSize;
f = (0:fftSize - 1)*deltaF;
f = mod(fftshift(f) + sampleRate/2, sampleRate) - sampleRate/2;

% Windows
hanning = hann(numPulseSamples);
hannWind = windowing('HanningModRisingDuration', numPulseSamples, sampleRate, 0.2);
rectWind = ones(numPulseSamples, 1);
hammWind = window(@hamming, numPulseSamples);

windows = [rectWind, hannWind, hanning, hammWind];

% FFT
fftWind = fftshift(fft(windows, fftSize), 1);

% Energy density
energyDensityWind = abs(fftWind).^2/sampleRate^2;

energyTime = sum(abs(windows).^2, 1)*deltaT;
energyFrec = sum(energyDensityWind, 1)*deltaF;

% Graphics: time and frequency domain
fig = figure;

ax1 = subplot(1, 2, 1);
plot(ax1, t, windows)

ax2 = subplot(1, 2, 2);
plot(ax2, f, energyDensityWind);
ax2.YLim = [0, 10];

% plot(ax2, f, 10*log10(energyDensityRect), f, 10*log10(energyDensityHann));
% ax2.YLim = [-30, 10];

ax2.XLim = [-4, 4];

path = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\';
name = 'window_example';
printfig(fig, path, name, {'emf', 'eps'});

% How condensed is the energy in the low frequencies
[fOrd, ind] = sort(abs(f));
cumEner = cumsum(energyDensityWind(ind, :))*deltaF;
% Normalize
cumEner = cumEner./repmat(energyTime, [fftSize, 1]);

ax = axes(figure);
plot(ax, fOrd, cumEner);
ax.XLim = [0, 10];


