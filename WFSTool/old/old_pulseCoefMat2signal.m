function [signal] = pulseCoefMat2signal(freq, coefMat, pulseLimits, sampleRate, startMarker, endMarker, type)
% coefMat. (numPulses x numChannels x numFreqs)
    
if nargin < 7
    type = 'sample';
end

startMarker = max(0, startMarker);
endMarker = min(pulseLimits(end), endMarker);

[newPulseLimits, indPulses] = cutPulseLimits(pulseLimits, startMarker, endMarker);

switch type
    case 'time' % Time version
        startSample = floor(startMarker*sampleRate);
        newPulseSampleLimits = floor(newPulseLimits*sampleRate);
    case 'sample' % Matlab sample version (indexing start at 1 in matlab)
        startSample = startMarker - 1;
        newPulseSampleLimits = newPulseLimits - 1;
end

coefMat = coefMat(indPulses, :, :);
signal = genSignal(coefMat, freq, newPulseSampleLimits-startSample, startSample, sampleRate);

end

function signal = genSignal(coefMat, freq, pulseLimits, offsetSample, sampleRate)
% coefMat. (numPulses x numChannels x numFreqs). The (i, j, k)-th element
% contains the complex coefficient of the i-th pulse for the j-th channel
% and the f-th frequency
% freq. numFreqs-element vector.
% pulseLimits. (numPulses x 2). The (i, 1)-th element is the sample where
% the i-th pulse starts. The (i, 2)-th element is the sample where the i-th
% pulse ends. The first sample of the signal is the sample 0, not 1.
% offsetSample. Scalar. The offset used to calculate the time vector.
% sampleRate. Scalar. Sample rate of the signal.

numPulses = size(coefMat, 1);
numChannels = size(coefMat, 2);
numFrequencies = size(coefMat, 3);

if isempty(pulseLimits)
    numSamples = 0;
else
    numSamples = pulseLimits(end);
end

signal = zeros(numSamples, numChannels);
for p = 1:numPulses
    startPulse = pulseLimits(p, 1);
    endPulse = pulseLimits(p, 2);
    numPulseSamples = endPulse - startPulse;
    
    ind = startPulse:endPulse-1;
    t = (ind + offsetSample)/sampleRate;
    
    for c = 1:numChannels
        aux = zeros(numPulseSamples, numFrequencies);
        for f = 1:numFrequencies
            amp = abs(coefMat(p, c, f));
            phase = angle(coefMat(p, c, f));
            aux(:, f) = amp*cos(2*pi*freq(f)*t + phase);
        end
        signal(ind + 1, c) = sum(aux, 2);
    end
end

end

