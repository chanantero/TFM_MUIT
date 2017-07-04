function [signal] = pulseCoefMat2signal(freq, coefMat, pulseLimits, sampleRate, startTime, endTime)
% coefMat. (numPulses x numChannels x numFreqs)

numPulses = size(coefMat, 1);
numChannels = size(coefMat, 2);
numFreqs = size(coefMat, 3);

startTime = max(0, startTime);
endTime = min(pulseLimits(end), endTime);

pulseLimits = reshape(pulseLimits', numPulses*2, 1);
indPulseLimits = find(pulseLimits >= startTime & pulseLimits <= endTime);
newPulseLimits = pulseLimits(indPulseLimits);

evenStart = rem(indPulseLimits(1), 2) == 0;
if evenStart
    newPulseLimits = [startTime; newPulseLimits];
end

oddEnd = rem(indPulseLimits(end), 2) == 1;
if oddEnd
    newPulseLimits = [newPulseLimits; endTime];
end

newPulseLimits = reshape(newPulseLimits, [2, numPulses])';
numNewPulses = size(newPulseLimits, 1);

startSample = floor(startTime*sampleRate);
endSample = floor(endTime*sampleRate);
newPulseSampleLimits = floor(newPulseLimits*sampleRate);
rel_pulseSampleLimits = newPulseSampleLimits - startSample;

numSamples = endSample - startSample + 1;

signal = zeros(numSamples, numChannels);
for p = 1:numNewPulses
    ind = (rel_pulseSampleLimits(p, 1):rel_pulseSampleLimits(p, 2)) + 1;
    t = (newPulseSampleLimits(p, 1):newPulseSampleLimits(p, 2))'/sampleRate;
    for c = 1:numChannels
        aux = zeros(numel(ind), numFreqs);
        for f = 1:numFreqs
            aux(:, f) = coefMat(p, c, f)*cos(2*pi*freq(f)*t);
        end
        signal(ind, c) = sum(aux, 2);
    end
end


end