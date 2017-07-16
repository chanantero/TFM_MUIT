function [signal, outOfRange] = pulseCoefMat2signal(freq, coefMat, pulseLimits, sampleRate, startMarker, endMarker, type)
% coefMat. (numPulses x numChannels x numFreqs)

if nargin < 7
    type = 'sample';
end

% startMarker = max(0, startMarker);
% endMarker = min(pulseLimits(end), endMarker);

switch type
    case 'time'
        startSample = floor(startMarker*sampleRate);
        endSample = floor(endMarker*sampleRate);
        pulseSampleLimits = floor(pulseLimits*sampleRate);
    case 'sample'
        startSample = startMarker - 1;
        endSample = endMarker - 1 + 1; % -1 because we work with 0 as the first index, not 1 as matlab does. +1 because we want the input of type "sample" to be the last sample to return
        pulseSampleLimits = pulseLimits - 1;
end

endSample = min(endSample, startMarker + 1e6); % Limit size of output

[cuttedPulseLimits, indPulses, outOfRange] = intervalSelection(pulseSampleLimits, startSample, endSample);

if isempty(indPulses)
    numChannels = size(coefMat, 2);
    signal = zeros(endSample - startSample, numChannels);
else
    startFirstPulse = pulseSampleLimits(indPulses(1),1);
    
    wholePulseLimits = pulseSampleLimits(indPulses, :) - startFirstPulse;
    
    coefMat = coefMat(indPulses, :, :);
    
    wholePulseSignal = genSignal(coefMat, freq, wholePulseLimits, startFirstPulse, sampleRate);
    
    % Map signal
    signal = zeros(endSample - startSample, size(wholePulseSignal, 2));
    selAbsInd = cuttedPulseLimits(1, 1):cuttedPulseLimits(end, 2) - 1;
    ind1 = selAbsInd - startSample + 1;
    ind2 = selAbsInd - startFirstPulse + 1;
    signal(ind1, :) = wholePulseSignal(ind2, :);
end



end

function signal = genSignal(coefMat, freq, pulseLimits, offsetSample, sampleRate)
% coefMat. (numPulses x numChannels x numFreqs). The (i, j, k)-th element
% contains the complex coefficient of the i-th pulse for the j-th channel
% and the f-th frequency
% freq. numFreqs-element vector.
% pulseLimits. (numPulses x 2) matrix.
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
    
    windMask = windowing2(numPulseSamples, sampleRate, 0.1);
    
    ind = startPulse:endPulse-1;
    t = (ind + offsetSample)/sampleRate;
    
    for c = 1:numChannels
        aux = zeros(numPulseSamples, numFrequencies);
        for f = 1:numFrequencies
            amp = abs(coefMat(p, c, f));
            phase = angle(coefMat(p, c, f));
            aux(:, f) = amp*cos(2*pi*freq(f)*t + phase);
        end
        signal(ind + 1, c) = sum(aux, 2).*windMask;
    end
end

end

function mask = windowing1(numSamples, constantRatio)
hanningRatio = 1 - constantRatio;
numHann = floor(hanningRatio*numSamples);
numPre = ceil(numHann/2);
numPost = numHann - numPre;
hannWind = hann(numHann);

preInd = 1:numPre;
postInd = (numSamples-numPost+1):numSamples;
constantInd = numPre+1:numSamples-numPost;

mask = zeros(numSamples, 1);
mask(preInd) = hannWind(1:numPre);
mask(postInd) = hannWind(numPre+1:end);
mask(constantInd) = 1;

end

function mask = windowing2(numSamples, sampleRate, risingDuration)
risingSamples = floor(sampleRate*risingDuration);

numPre = min(risingSamples, floor(numSamples/2));
numPost = numPre;
hannWind = hann(numPre*2);

preInd = 1:numPre;
postInd = (numSamples-numPost+1):numSamples;
constantInd = numPre+1:numSamples-numPost;

mask = zeros(numSamples, 1);
mask(preInd) = hannWind(1:numPre);
mask(postInd) = hannWind(numPre+1:end);
mask(constantInd) = 1;

end

