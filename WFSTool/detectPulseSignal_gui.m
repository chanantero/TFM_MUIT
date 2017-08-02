function [coeff, correspondingIndexes, yPulseLimits] = detectPulseSignal_gui(ax, xPulseLimits, y, ySampleRate, marginRatio)
% xPulseLimits. (numPulses x 2) matrix. The first sample has the index 0.
% y. (numSamples x numChannels) matrix
% ySampleRate. Scalar. Natural number
% marginRatio. Scalar between 0 and 1

numPulses = size(xPulseLimits, 1);

% Correlate a mask created with xPulseLimits with the received signal y

% Generate the basic mask
xPulseLimitsSamp = floor(xPulseLimits * ySampleRate);
% beggining = xPulseLimitsSamp(1, 1);
ending = xPulseLimitsSamp(end, 2);
numSampMask = ending;
mask = zeros(numSampMask, 1);
for p = 1:numPulses
    mask(xPulseLimitsSamp(p, 1)+1:xPulseLimitsSamp(p, 2)) = 1;
end

% Correlate it
meanDur = mean(diff(xPulseLimits, 1, 2));
step = floor(meanDur * 0.1 * ySampleRate);
% lowLimit = -(numSampMask - 1);
% highLimit = numSamples-1;
% Estimate a range of shift positions of the mask. The range will be guided
% by a concrete number of pulses. We assume that the shift of the received
% signal y is less than a concrete number of pulses.
NP = 10;
lowLimit = -xPulseLimitsSamp(NP, 2);
highLimit = xPulseLimitsSamp(NP, 2);

% Apply mask only to a fragment of the signal
y_ = y(1:min(highLimit*2, end), :);
numSamples_ = size(y_, 1);

firstIteration = true;
while step > 1
    shiftPos = highLimit:-step:lowLimit;
    N = numel(shiftPos);
    corr = zeros(N, 1); % Correlation
    for k = 1:N
        % Shift mask
        shiftedMask = shiftAndCrop(mask, shiftPos(k), numSamples_);
        
        % Apply mask
        corr(k, :) = sum(abs(y_).*shiftedMask);
        
%         fprintf('%d/%d\n', k, N)
    end
    
    % Normalize the correlation values
    corr = corr/max(corr);
    
    if firstIteration
        % Find the peaks in the correlation
        delta = (max(corr) - min(corr)) * 0.1;
        [maxtab, ~] = peakdet(corr, delta);
        
        % The maximum peak has the value 1. Find all the peaks that are above
        % 0.985 and select the first one
        cand = find(maxtab(:, 2) >= 0.985, 1, 'first');
        ind = maxtab(cand, 1);
        
        firstIteration = false;
    else
        % Just find the maximum value
        [~, ind] = max(corr);
    end
    
    % Set the new step
    if ind == lowLimit || ind == highLimit
        break;
    end
    
    lowLimit = shiftPos(ind + 1);
    highLimit = shiftPos(ind - 1);
    step = floor((highLimit - lowLimit)*0.1);
end

shiftSamp = shiftPos(ind);

% Find the pulse signal coefficients
% The new pulse limits are the original ones shifted
yPulseLimitsSamp = xPulseLimitsSamp + shiftSamp;
valid = find(all(yPulseLimitsSamp >= 0, 2));
numDetPulses = numel(valid);

yPulseLimits = yPulseLimitsSamp(valid, :);
correspondingIndexes = valid;
coeff = zeros(numDetPulses, 1);
for k = 1:numDetPulses
    % Discard the extremes of the pulse
    startPulse = yPulseLimitsSamp(valid(k), 1);
    endPulse = yPulseLimitsSamp(valid(k), 2);
    durPulse = endPulse - startPulse;
    margin = floor(durPulse * marginRatio);
    ind = (startPulse + margin):(endPulse - 1 - margin);
    coeff(k) = mean(y(ind));
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