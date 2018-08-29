function [varargout] = getFrequencyFilter( magnitudeFilterOrder, phaseFilterOrder, Fs, varargin )

magnitudeFlag = ~isempty(magnitudeFilterOrder);
phaseFlag = ~isempty(phaseFilterOrder);

p = inputParser;
addParameter(p, 'analytical', false)
parse(p, varargin{:})
analyticalFlag = p.Results.analytical;

if analyticalFlag
    if magnitudeFlag
        [hFreqDepend, delayMagnFilter] = magnitudeFilterAnalytical(magnitudeFilterOrder, Fs);
    end
    if phaseFlag
        [hPhase, delayPhaseFilter] = phaseFilterAnalytical(phaseFilterOrder, Fs);
    end
    
    if magnitudeFlag && phaseFlag
        delay = delayMagnFilter + delayPhaseFilter;
        hTotal = conv(hFreqDepend, hPhase);
        varargout = {hTotal, delay, hFreqDepend, delayMagnFilter, hPhase, delayPhaseFilter};
    else
        if magnitudeFlag
            varargout = {hFreqDepend, delayMagnFilter};
        elseif phaseFlag
            varargout = {hPhase, delayPhaseFilter};
        else
            varargout = {};
        end
    end
else
    if magnitudeFlag
        % Filter 1: Filter for the frequency dependence: sqrt(k)
        [hFreqDepend, delayMagnFilter] = magnitudeFilter(magnitudeFilterOrder, Fs);
    end
    
    if phaseFlag
        % Filter 2: Hilbert filter
        [hHilbert, delayHilbertFilter] = HilbertFilter(phaseFilterOrder, Fs);
    end
    
    if magnitudeFlag && phaseFlag
        % Total filter to apply to the signal.
        delta = zeros(size(hHilbert));
        delta(delayHilbertFilter + 1) = 1;
        hTotal = 1/sqrt(2)*conv((delta - hHilbert), hFreqDepend);
        
        delay = delayMagnFilter + delayHilbertFilter;
        
        varargout = {hTotal, delay, hFreqDepend, delayMagnFilter, hHilbert, delayHilbertFilter};
    else
        if magnitudeFlag
            varargout = {hFreqDepend, delayMagnFilter};
        elseif phaseFlag
            varargout = {hHilbert, delayHilbertFilter};
        else
            varargout = {};
        end
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

function [hMag, delay] = magnitudeFilterAnalytical(magnitudeFilterOrder, Fs)
% Frank Schutz. (2015). Sound Field Synthesis for Line Source Array Applications in Large-Scale Sound Reinforcement. Dissertation.
if rem(magnitudeFilterOrder, 2) == 0
    magnitudeFilterOrder = magnitudeFilterOrder + 1;
end

n = -(magnitudeFilterOrder-1)/2:(magnitudeFilterOrder-1)/2;
hMag = -fresnels(sqrt(2*n))./(sqrt(2*pi)*n.^(3/2));
zeroInd = (magnitudeFilterOrder-1)/2 + 1;
hMag(zeroInd) = 2/3*sqrt(pi);

c = 340;
hMag = hMag*2/5*sqrt(Fs/c);

delay = (zeroInd - 1);
end

function [hPha, delay] = phaseFilterAnalytical(phaseFilterOrder, Fs)
% Frank Schutz. (2015). Sound Field Synthesis for Line Source Array Applications in Large-Scale Sound Reinforcement. Dissertation.

if rem(phaseFilterOrder, 2) == 0
    phaseFilterOrder = phaseFilterOrder + 1;
end

n = -(phaseFilterOrder-1)/2:(phaseFilterOrder-1)/2;
nEvenInd = rem(n,2) == 0;
nOddInd = rem(n,2) ~= 0;
zeroInd = find(n == 0);
hPha = zeros(1, phaseFilterOrder);
hPha(nEvenInd) = 0;
hPha(nOddInd) = -sqrt(2)./(pi*n(nOddInd));
hPha(zeroInd) = 1/sqrt(2);

delay = (zeroInd - 1);

end
