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
    % Form D)
    % Es eficiente, permite el solapamiento y genera bien la máscara
        % Generate signal
    startFirstPulse = cuttedPulseLimits(1, 1);
    coefMat = coefMat(indPulses, :, :);
    offsetCuttedPulseLimits = cuttedPulseLimits - startFirstPulse;
    pulseSignal = genSignalNew(coefMat, freq, offsetCuttedPulseLimits, pulseSampleLimits(indPulses, :), startFirstPulse, sampleRate);
        % Map signal
    signal = zeros(endSample - startSample, size(pulseSignal, 2));
    selAbsInd = cuttedPulseLimits(1, 1):cuttedPulseLimits(end, 2) - 1;
    ind = selAbsInd - startSample + 1;
    signal(ind, :) = pulseSignal;
    
    
%     % Form C)
%     % Es eficiente, permite el solapamiento, pero no genera bien la máscara
%         % Generate signal
%     startFirstPulse = cuttedPulseLimits(1, 1);
%     coefMat = coefMat(indPulses, :, :);
%     offsetCuttedPulseLimits = cuttedPulseLimits - startFirstPulse;
%     pulseSignal = genSignal(coefMat, freq, offsetCuttedPulseLimits, startFirstPulse, sampleRate);
%             % Map signal
%     signal = zeros(endSample - startSample, size(pulseSignal, 2));
%     selAbsInd = cuttedPulseLimits(1, 1):cuttedPulseLimits(end, 2) - 1;
%     ind = selAbsInd - startSample + 1;
%     signal1(ind, :) = pulseSignal;
    
%     % Form B)
%     % Es eficiente, genera bien la máscara, pero no permite solapamiento
%         % Generate signal
%     startFirstPulse = cuttedPulseLimits(1, 1);
%     coefMat = coefMat(indPulses, :, :);
%     offsetCuttedPulseLimits = cuttedPulseLimits - startFirstPulse;
%     pulseSignal = genSignalNoMask(coefMat, freq, offsetCuttedPulseLimits, startFirstPulse, sampleRate);
%             % Map signal
%     signal = zeros(endSample - startSample, size(pulseSignal, 2));
%     selAbsInd = cuttedPulseLimits(1, 1):cuttedPulseLimits(end, 2) - 1;
%     ind = selAbsInd - startSample + 1;
%     signal(ind, :) = pulseSignal;
%     
%         % Generate mask
%     risingDuration = 0.1;
%     wholePulseSampleLimits = pulseSampleLimits(indPulses, :);
%     mask = genMask(wholePulseSampleLimits - wholePulseSampleLimits(1,1), sampleRate, risingDuration);    % Map masking
%     masking = zeros(endSample - startSample, 1);
%     indA = selAbsInd - startSample + 1;
%     indB = selAbsInd - wholePulseSampleLimits(1,1) + 1;
%     masking(indA) = mask(indB);
% 
%         % Apply mask to signal
%     numChann = size(signal, 2);
%     signal2 = signal.*repmat(masking, [1 numChann]);
    
    
%     % Form A)
%     % Nada eficiente, pero genera bien la máscara y permite solapamiento
%     startFirstPulse = pulseSampleLimits(indPulses(1),1);
%     
%     wholePulseLimits = pulseSampleLimits(indPulses, :) - startFirstPulse;
%     
%     coefMat = coefMat(indPulses, :, :);
%     
%     tic
%     wholePulseSignal = genSignal(coefMat, freq, wholePulseLimits, startFirstPulse, sampleRate);
%     toc
%     
%     % Map signal
%     signal = zeros(endSample - startSample, size(wholePulseSignal, 2));
%     selAbsInd = cuttedPulseLimits(1, 1):cuttedPulseLimits(end, 2) - 1;
%     ind1 = selAbsInd - startSample + 1;
%     ind2 = selAbsInd - startFirstPulse + 1;
%     signal(ind1, :) = wholePulseSignal(ind2, :);
end

end

function mask = genMask(pulseSampleLimits, sampleRate, risingDuration)
    
mask = zeros(pulseSampleLimits(end) - 1, 1);
numPulses = size(pulseSampleLimits, 1);

for p = 1:numPulses
    startSample = pulseSampleLimits(p, 1);
    endSample = pulseSampleLimits(p, 2);
    numPulseSamples = endSample - startSample;
    windMask = windowing('HanningModRisingDuration', numPulseSamples, sampleRate, risingDuration);
    
    ind = (startSample:endSample-1) + 1;
    mask(ind) = windMask;
end

end

function signal = genSignalNew(coefMat, freq, pulseLimits, wholePulseLimits, offsetSample, sampleRate)
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
    
    ind = startPulse:endPulse-1;
    t = (ind + offsetSample)'/sampleRate;
    
    %     windMask = ones(numPulseSamples, 1); % Uniform window
    indSampMask = wholePulseLimits(p,1):wholePulseLimits(p,2)-1;
    windMask = windowing('HanningModRisingDuration', length(indSampMask), sampleRate, 0.1);
    [flag, flagInd] = ismember(indSampMask, ind + offsetSample);
    windMaskAdapted = zeros(size(t));
    windMaskAdapted(flagInd(flag)) = windMask(flag);
    
    iszero = coefMat(p, :, :) == 0;
    nonZeroChannels = find(any(~iszero, 3));
    numNonZeroChannels = length(nonZeroChannels);
    for c = 1:numNonZeroChannels
        pulse = genTones(coefMat(p, nonZeroChannels(c), :), freq, t);
        signal(ind + 1, nonZeroChannels(c)) = signal(ind + 1, nonZeroChannels(c)) + pulse.*windMaskAdapted; % Si solo está pulse.*windMask, no se pueden solapar los pulsos. Con signal(ind + 1, c) + pulse.*windMask, sí puede haber solapamiento.
    end
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

function signal = genSignalNoMask(coefMat, freq, pulseLimits, offsetSample, sampleRate)
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
    
%     numPulseSamples = endPulse - startPulse;
%     windMask = ones(numPulseSamples, 1); % Uniform window
%     windMask = windowing('HanningModRisingDuration', numPulseSamples, sampleRate, 0.1);
    
    ind = startPulse:endPulse-1;
    t = (ind + offsetSample)/sampleRate;
    
    iszero = coefMat(p, :, :) == 0;
    nonZeroChannels = find(any(~iszero, 3));
    numNonZeroChannels = length(nonZeroChannels);
    for c = 1:numNonZeroChannels
        pulse = genTones(coefMat(p, nonZeroChannels(c), :), freq, t);
        
%         signal(ind + 1, nonZeroChannels(c)) = signal(ind + 1, nonZeroChannels(c)) + pulse.*windMask; % Si solo está pulse.*windMask, no se pueden solapar los pulsos. Con signal(ind + 1, c) + pulse.*windMask, sí puede haber solapamiento.
        signal(ind + 1, nonZeroChannels(c)) = signal(ind + 1, nonZeroChannels(c)) + pulse;
    end
end

end

function pulse = genTones(coefficients, frequencies, time)

nonzero = find(coefficients ~= 0);
frequencies = frequencies(nonzero);
coefficients = coefficients(nonzero);

numFrequencies = numel(frequencies);
time = time(:);
numSamples = length(time);

memoryFriendly = false;
try
    aux = zeros(numSamples, numFrequencies);
catch ME
    memoryFriendly = true;
end

if memoryFriendly
    pulse = zeros(numSamples, 1);
    for f = 1:numFrequencies
            amp = abs(coefficients(f));
            phase = angle(coefficients(f));
            pulse = pulse + amp*cos(2*pi*frequencies(f)*time + phase);
    end
else
    for f = 1:numFrequencies
            amp = abs(coefficients(f));
            phase = angle(coefficients(f));
            aux(:, f) = amp*cos(2*pi*frequencies(f)*time + phase);
    end
    pulse = sum(aux, 2);
end



end




