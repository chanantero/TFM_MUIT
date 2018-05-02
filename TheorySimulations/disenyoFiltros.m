%% Creaci�n de respuestas al impulso para el filtrado temporal de se�ales

% Filter 1: Filter for the frequency dependence: sqrt(k)
    % Parameters
c = 340; % Sound speed
loudspeakerSeparation = 0.18;
Fs = 44100; % Sampling frequency 
fmaxAliasing = c/(2*loudspeakerSeparation); % Maximum frequency free from aliasing

    % Filter desing
maxSyntFreq = fmaxAliasing*2; % Maximum frequency at which we are interested in synthesizing. Above that frequency the filter response is irrelevant, since there will be aliasing and hence, we won't use signals with such high frequencies
fNyquist = Fs/2;
freqs = 0:40:maxSyntFreq;
% Modify the vector of frequencies since the function firls requires that
freqs_aux = [freqs(1); reshape([freqs(2:end-1); freqs(2:end-1)], [(numel(freqs) - 2)*2, 1]); freqs(end)];
magnFilter = sqrt(1i*freqs_aux/c);
freqs_aux = [freqs_aux; freqs_aux(end); fNyquist];
magnFilter = [magnFilter; magnFilter(end); magnFilter(end)];
% plot(freqs_aux, abs(magnFilter))
ord=512;
hFreqDepend = firls(ord, freqs_aux/fNyquist, abs(magnFilter));
fvtool(hFreqDepend, 1, 'OverlayedAnalysis', 'phase')

[~, indMax] = max(hFreqDepend);
delayFilter1 = indMax - 1;

% Filter 2: Hilbert filter
N = 16384;
d = fdesign.hilbert('N,TW', N, 20, Fs);
hd = design(d, 'equiripple', 'SystemObject', true);
% zerophase(hd,-fs/2:0.1:fs/2,fs)
hHilbert = hd.Numerator;
fvtool(hHilbert, 1, 'OverlayedAnalysis', 'phase')

delayFilter2 = N/2;

% Total filter to apply to the signal.
% Notice that this filter needs to be applied only once, to the signal from
% the source. That gives an advantage in the computing complexity
delta = zeros(size(hHilbert));
delta(delayFilter2) = 1;
hTotal = 1/sqrt(2)*conv((delta - hHilbert), hFreqDepend);
fvtool(hTotal, 1, 'OverlayedAnalysis', 'phase')

save('WFSfilter.mat', 'hTotal', 'Fs')

f = 900;
dt = 1/Fs;
t = 0:dt:2;
x = cos(2*pi*f*t);
y1 = filter(hFreqDepend, 1, x);
y2 = filter(hHilbert, 1, x);
yTotal1 = filter(hHilbert, 1, y1);
yTotal = filter(hTotal, 1, x);

y1_undelayed = [y1(delayFilter1+1:end), zeros(1, delayFilter1)]; 
plot(t, x, t, y1_undelayed)

y2_undelayed = [y2(delayFilter2+1:end), zeros(1, delayFilter2)]; 
plot(t, x, t, y2_undelayed)

yTotal1_undelayed = [yTotal1(delayFilter1 + delayFilter2 + 1:end), zeros(1, delayFilter1 + delayFilter2)];
plot(t, x, t, yTotal1_undelayed)

yTotal_undelayed = [yTotal(delayFilter1 + delayFilter2 + 1:end), zeros(1, delayFilter1 + delayFilter2)];
plot(t, x, t, yTotal_undelayed)


plot(t, x, t, y1)
plot(t, x, t, y2)

%% Various frequencies
f = [900 440];
dt = 1/Fs;
t = (0:dt:2)';
x = cos(2*pi*repmat(f, [numel(t), 1]).*repmat(t, [1, numel(f)]));
yTotal = filter(hTotal, 1, x);

yTotal_undelayed = [yTotal(delayFilter1 + delayFilter2 + 1:end, :); zeros(delayFilter1 + delayFilter2, numel(f))];
plot(t, x, t, yTotal_undelayed)

%% Tests
b = zeros(1, 10); b(4) = 1;
x = ones(1, 15);
zi = ones(1, 9);
y = filter(b, 1, x, zi);

ax = axes(figure);
plot(ax, y)

N = length(hTotal);
Htotal = fft(hTotal);
df = 1/N;
f = (0:df:1-df)*44100;
f_fftshift = (-floor(N/2):floor(N/2))*44100/N;
ax = axes(figure, 'NextPlot', 'Add');
plot(ax, f_fftshift, abs(fftshift(Htotal)))
ax.YLim = [0, 4];
ax.XLim = [-2500, 2500];
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = '|H(f)|';
c = 340;
plot(ax, f_fftshift, sqrt(abs(f_fftshift)/343))
legend(ax, 'Real', 'Ideal')
ax.Title.String = 'Frequency response';

t = (0:N-1)/Fs;
ax = axes(figure);
plot(ax, t, abs(hTotal))
ax.XLim = [0, max(t)];
ax.XLabel.String = 'Time (s)';
ax.YLabel.String = '|h(t)|';
ax.Title.String = 'Impulse response';
ax.YLim = [-0.1, 3];

N = length(hFreqDepend);
HfreqDepend = fft(hFreqDepend);
df = 1/N;
f = (0:df:1-df)*44100;
f_fftshift = (-floor(N/2):floor(N/2))*44100/N;
ax = axes(figure, 'NextPlot', 'Add');
plot(ax, f_fftshift, abs(fftshift(HfreqDepend)))
ax.YLim = [0, 4];
ax.XLim = [-2500, 2500];
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = '|H(f)|';
c = 340;
plot(ax, f_fftshift, sqrt(f_fftshift/343))
plot(ax, freqs_aux, abs(magnFilter))

N = length(hHilbert);
Hhilbert = fft(hHilbert);
df = 1/N;
f = (0:df:1-df)*44100;
f_fftshift = (-floor(N/2):floor(N/2))*44100/N;
ax = axes(figure, 'NextPlot', 'Add');

plot(ax, f_fftshift, abs(fftshift(Hhilbert)))
ax.YLim = [0, 4];
ax.XLim = [-2500, 2500];

yyaxis right
phase = unwrap(angle(fftshift(Hhilbert)));
plot(ax, f_fftshift, phase)
ax.YLim = [-40000, 40000];
ax.XLim = [-2500, 2500];

ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = '|H(f)|';

