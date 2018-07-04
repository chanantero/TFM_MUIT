function [varargout] = getFrequencyFilter( magnitudeFilterOrder, hilbertFilterOrder, Fs )

magnitudeFlag = ~isempty(magnitudeFilterOrder);
hilbertFlag = ~isempty(hilbertFilterOrder);

if magnitudeFlag
    % Filter 1: Filter for the frequency dependence: sqrt(k)
    [hFreqDepend, delayMagnFilter] = magnitudeFilter(magnitudeFilterOrder, Fs);
end

if hilbertFlag
    % Filter 2: Hilbert filter
    [hHilbert, delayHilbertFilter] = HilbertFilter(hilbertFilterOrder, Fs);
end

if magnitudeFlag && hilbertFlag
    % Total filter to apply to the signal.
    delta = zeros(size(hHilbert));
    delta(delayHilbertFilter + 1) = 1;
    hTotal = 1/sqrt(2)*conv((delta - hHilbert), hFreqDepend);
    
    delay = delayMagnFilter + delayHilbertFilter;
    
    varargout = {hTotal, delay, hFreqDepend, delayMagnFilter, hHilbert, delayHilbertFilter};
else
    if magnitudeFlag
        varargout = {hFreqDepend, delayMagnFilter};
    elseif hilbertFlag
        varargout = {hHilbert, delayHilbertFilter};
    else
        varargout = {};
    end
end

end

function [hFreqDepend, delay] = magnitudeFilter(magnitudeFilterOrder, Fs)
% Parameters
c = 340; % Sound speed
loudspeakerSeparation = 0.18;
% Fs = 44100; % Sampling frequency
fmaxAliasing = c/(2*loudspeakerSeparation); % Maximum frequency free from aliasing

% Filter desing
maxSyntFreq = fmaxAliasing*1.1; % Maximum frequency at which we are interested in synthesizing. Above that frequency the filter response is irrelevant, since there will be aliasing and hence, we won't use signals with such high frequencies
fNyquist = Fs/2;
freqs = 0:40:maxSyntFreq;
% Modify the vector of frequencies since the function firls requires that
freqs_aux = [freqs(1); reshape([freqs(2:end-1); freqs(2:end-1)], [(numel(freqs) - 2)*2, 1]); freqs(end)];
magnFilter = sqrt(freqs_aux/c); % abs(sqrt(1i*freqs_aux/c))
freqs_aux = [freqs_aux; freqs_aux(end); fNyquist];
magnFilter = [magnFilter; magnFilter(end); magnFilter(end)];
hFreqDepend = firls(magnitudeFilterOrder, freqs_aux/fNyquist, magnFilter);

[~, indMax] = max(hFreqDepend);
delay = indMax - 1;
end

function [hHilbert, delay] = HilbertFilter(hilbertFilterOrder, Fs)
d = fdesign.hilbert('N,TW', hilbertFilterOrder, 20, Fs);
hd = design(d, 'equiripple', 'SystemObject', true);
hHilbert = hd.Numerator;

delay = hilbertFilterOrder/2;
end