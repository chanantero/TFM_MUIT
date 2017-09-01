% Diapositiva Generación Señal Pulsos

% Parameters
pulseCoefMat = ...
    [1i;
    0;
    0;
    0];

pulseLimits = ...
    [0 2;
    3 4;
    5 7
    9 10];

freqs = 2;
sampleRate = 10000;

% Apply window. Generate signal
signal = pulseCoefMat2signal(pulseCoefMat, pulseLimits, freqs, sampleRate, 0, 10, 'type_marker', 'time', 'type_pulseLimits', 'time');
numSamples = size(signal, 1);
t = (0:numSamples-1)/sampleRate;

% Base Tone
tone = cos(2*pi*t*freqs(1));

% Apply coefficient
toneCoef = real(pulseCoefMat(1,1,1)*exp(1i*2*pi*freqs(1)*t));

% Multiple tones. 1 frequency, 1 channel
pulseCoefMat = ...
    [1i;
    0;
    2
    0];

freqs = 2;

multPulses = pulseCoefMat2signal(pulseCoefMat, pulseLimits, freqs, sampleRate, 0, 10, 'type_marker', 'time', 'type_pulseLimits', 'time');

% Multiple tones. 2 frequencies, 1 channel
pulseCoefMat = ...
    [1i;
    0;
    2;
    0];

pulseCoefMat(:, 1, 2) = ...
    [-1i; 0; 0; 0];

freqs = [2, 0.25];

multFreq = pulseCoefMat2signal(pulseCoefMat, pulseLimits, freqs, sampleRate, 0, 10, 'type_marker', 'time', 'type_pulseLimits', 'time');

% Multiple channels
pulseCoefMat = ...
    [1i 0;
    0 1;
    2 0;
    0 -0.5];

pulseCoefMat(:, :, 2) = ...
    [-1i 0; 0 0; 0 1; 0 0];

freqs = [2, 0.25];

complete = pulseCoefMat2signal(pulseCoefMat, pulseLimits, freqs, sampleRate, 0, 10, 'type_marker', 'time', 'type_pulseLimits', 'time');


% Graphs
path = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\';

fig = figure;
ax = axes(fig);
ax.YLim = [-2, 2];
ax.XLabel.String = 'Time (s)';
l = line(ax, t, zeros(size(t)));

l.YData = tone;
printfig(fig, path, 'baseTone', {'emf'})

l.YData = toneCoef;
printfig(fig, path, 'coefTone', {'emf'})

l.YData = signal;
printfig(fig, path, 'singlePulse', {'emf'})

l.YData = multPulses;
printfig(fig, path, 'multPulses', {'emf'})

l.YData = multFreq;
printfig(fig, path, 'multFreq', {'emf'})

l = line(ax, t, complete);
printfig(fig, path, 'completePulseSignal', {'emf'})

