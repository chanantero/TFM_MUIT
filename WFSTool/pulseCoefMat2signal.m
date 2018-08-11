function [signal, outOfRange] = pulseCoefMat2signal(coefMat, pulseLimits, freq, sampleRate, varargin)
% coefMat. (numPulses x numChannels x numFreqs)

p = inputParser;

addOptional(p, 'startMarker', [])
addOptional(p, 'endMarker', [])
addParameter(p, 'type_pulseLimits', 'sample')
addParameter(p, 'type_marker', 'sample')

parse(p, varargin{:});

startMarker = p.Results.startMarker;
endMarker = p.Results.endMarker;
type_pulseLimits = p.Results.type_pulseLimits;
type_marker = p.Results.type_marker;

switch type_pulseLimits
    case 'time'
        pulseSampleLimits = floor(pulseLimits*sampleRate);
    case 'sample'
        pulseSampleLimits = pulseLimits;
end

if ismember('startMarker', p.UsingDefaults) % Assume that endMarker wasn't introduced either
    type_marker = 'sample';
    startMarker = min(pulseSampleLimits(:)) + 1;
    endMarker = max(pulseSampleLimits(:)) + 1 - 1; % +1 because endMarker must be a sample index starting at 1, not as 0 as with Matlab. -1 because we want endMarker to be the last sample to return
end

% startMarker = max(0, startMarker);
% endMarker = min(pulseLimits(end), endMarker);

switch type_marker
    case 'time'
        startSample = floor(startMarker*sampleRate);
        endSample = floor(endMarker*sampleRate);
    case 'sample'
        startSample = startMarker - 1;
        endSample = endMarker - 1 + 1; % -1 because we work with 0 as the first index, not 1 as matlab does. +1 because we want the input of type "sample" to be the last sample to return
end

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
    
%     windMask = ones(numPulseSamples, 1); % Uniform window
    windMask = windowing('HanningModRisingDuration', numPulseSamples, sampleRate, 0.1);
    
    ind = startPulse:endPulse-1;
    t = (ind + offsetSample)/sampleRate;
    
    iszero = coefMat(p, :, :) == 0;
    nonZeroChannels = find(any(~iszero, 3));
    numNonZeroChannels = length(nonZeroChannels);
    for c = 1:numNonZeroChannels
        pulse = genTones(coefMat(p, nonZeroChannels(c), :), freq, t);
        
        signal(ind + 1, nonZeroChannels(c)) = signal(ind + 1, nonZeroChannels(c)) + pulse.*windMask; % Si solo está pulse.*windMask, no se pueden solapar los pulsos. Con signal(ind + 1, c) + pulse.*windMask, sí puede haber solapamiento.
    end
end

end

function pulse = genTones(coefficients, frequencies, time)
numFrequencies = numel(frequencies);
numSamples = numel(time);
time = time(:);

memoryFriendly = false;
try
    aux = zeros(numSamples, numFrequencies);
catch ME
    memoryFriendly = true;
end

iszero = coefficients == 0;

if memoryFriendly
    pulse = zeros(numSamples, 1);
    for f = 1:numFrequencies
        if ~iszero(f)
            amp = abs(coefficients(f));
            phase = angle(coefficients(f));
            pulse = pulse + amp*cos(2*pi*frequencies(f)*time + phase);
        end
    end
else
    for f = 1:numFrequencies
        if ~iszero(f)
            amp = abs(coefficients(f));
            phase = angle(coefficients(f));
            aux(:, f) = amp*cos(2*pi*frequencies(f)*time + phase);
        end
    end
    pulse = sum(aux, 2);
end



end




