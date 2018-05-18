function [hTotal, delay] = getFrequencyFilter( magnitudeFilterOrder, hilbertFilterOrder, Fs )

% Filter 1: Filter for the frequency dependence: sqrt(k)
    % Parameters
c = 340; % Sound speed
loudspeakerSeparation = 0.18;
% Fs = 44100; % Sampling frequency 
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
hFreqDepend = firls(magnitudeFilterOrder, freqs_aux/fNyquist, abs(magnFilter));

[~, indMax] = max(hFreqDepend);
delayFilter1 = indMax - 1;

% Filter 2: Hilbert filter
d = fdesign.hilbert('N,TW', hilbertFilterOrder, 20, Fs);
hd = design(d, 'equiripple', 'SystemObject', true);
hHilbert = hd.Numerator;

delayFilter2 = hilbertFilterOrder/2;

% Total filter to apply to the signal.
delta = zeros(size(hHilbert));
delta(delayFilter2) = 1;
hTotal = 1/sqrt(2)*conv((delta - hHilbert), hFreqDepend);

delay = delayFilter1 + delayFilter2 - 1;

end

