function [coeff, correspondingIndexes] = detectPulseSignal(xPulseLimits, y, ySampleRate, marginRatio)
% xPulseLimits. (numPulses x 2) matrix
% y. (numSamples x numChannels) matrix
% ySampleRate. Scalar. Natural number
% marginRatio. Scalar between 0 and 1

numPulses = size(xPulseLimits, 1);
numSamples = size(y, 1);
numChannels = size(y, 2);

% Correlate a mask created with xPulseLimits with the received signal y

% Generate the basic mask
xPulseLimitsSamp = floor(xPulseLimits * ySampleRate);
beggining = xPulseLimitsSamp(1, 1);
ending = xPulseLimitsSamp(end, 2);
numSampMask = beggining - ending + 1;
mask = zeros(numSampMask, 1);
for p = 1:numPulses
    mask(xPulseLimitsSamp(p, 1)+1:xPulseLimitsSamp(p, 2)) = 1;
end

% Correlate it
N = numSamples + numSampMask - 1;
shiftPos = numSamples - 1:-1:-(numSampMask - 1);
corr = zeros(N, numChannels); % Correlation
for k = 1:N
    % Shift mask
    shiftedMask = shiftAndCrop(mask, shiftPos(k), numSamples);
    
    % Apply mask
    corr(k, :) = sum(y.*repmat(shiftedMask, 1, numChannels), 1);
end

% Find the indexes of maximum correlation
[~, maxInd] = max(corr, [], 1);

% Find the pulse signal coefficients
coeff = cell(numChannels, 1);
correspondingIndexes = cell(numChannels, 1);
for cy = 1:numChannels
    % The new pulse limits are the original ones shifted
    yPulseLimitsSamp = xPulseLimitsSamp + shiftPos(maxInd(cy));
    valid = find(all(yPulseLimitsSamp >= 0, 2));        
    numDetPulses = numel(valid);
    
    correspondingIndexes{cy} = valid;
    coeff{cy} = zeros(numDetPulses, 1);
    for k = 1:numDetPulses
        % Discard the extremes of the pulse
        startPulse = yPulseLimitsSamp(valid(k), 1);
        endPulse = yPulseLimitsSamp(valid(k), 2);
        durPulse = endPulse - startPulse;
        margin = floor(durPulse * marginRatio);
        ind = (startPulse + margin):(endPulse - 1 - margin);
        coeff{cy}(k) = mean(x(ind));
    end
    
end

end

function y = shiftAndCrop(x, shiftPos, newLength)

N = size(x, 1);

y = zeros(newLength, 1);

firstYind = max(1, 1 + shiftPos);
lastYind = min(newLength, N + shiftPos);

firstXind = max(1, 1 - shiftPos);
lastXind = min(N, newLength - shiftPos);

y(firstYind:lastYind) = x(firstXind:lastXind);

end