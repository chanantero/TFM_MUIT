function [coeff, indicesFound, yPulseLimits] = detectPulseSignal(xPulseLimits, y, ySampleRate, marginType, marginValue)
% xPulseLimits. The first sample has the index 0.
% y. (numSamples x numChannels) matrix
% ySampleRate. Scalar. Natural number
% marginRatio. Scalar between 0 and 1

% Parse inputs
if nargin == 3
    marginType = 'marginRatio';
    marginValue = 0.1;
end    

% Calculate basic parameters
numPulses = size(xPulseLimits, 1);
numSamples = size(y, 1);

% Find the temporal shift.
% Correlate a mask created with xPulseLimits with the received signal y and
% find the maximum value of the correlation

    % Generate the mask limits
    xPulseLimitsSamp = floor(xPulseLimits * ySampleRate);
    startingSamples = xPulseLimitsSamp(:, 1) + 1; % +1 because xPulseLimitsSamp have 0 as the first sample
    endingSamples = xPulseLimitsSamp(:, 2) + 1 -1; % +1 because xPulseLimitsSamp have 0 as the first sample. -1 because the ending indices indicate an open interval (they are not included in the intervals)
    maskLimits = [startingSamples, endingSamples];

    % Generate the minimum and maximum lag and the resolution step for the
    % first iteration of the correlation
    meanDur = mean(diff(xPulseLimits, 1, 2));
    firstStep = floor(meanDur * 0.1 * ySampleRate);
    % Estimate a range of shift positions of the mask. The range will be guided
    % by a concrete number of pulses. We assume that the shift of the received
    % signal y is less than a concrete number of pulses.
    NP = min(10, numPulses);
    minLag = -xPulseLimitsSamp(NP, 2);
    maxLag = max(xPulseLimitsSamp(NP, 2), floor(numel(y)/2));

    % Correlate to find the shift
    shiftSamp = findShift(abs(y), maskLimits, firstStep, minLag, maxLag);


% Find the pulse signal coefficients
    
    % Shift the pulse limits
    yPulseLimitsSamp = xPulseLimitsSamp + shiftSamp; % The new pulse limits are the original ones shifted
    indicesFound = find(all(yPulseLimitsSamp >= 0 & yPulseLimitsSamp <= numSamples, 2)); % Don't consider the pulses that are not complete, i.e., they are out of the available signal
    numDetPulses = numel(indicesFound);
    yPulseLimits = yPulseLimitsSamp(indicesFound, :);
    startPulses = yPulseLimits(:, 1);
    endPulses = yPulseLimits(:, 2);

    % Margin. Calculate the number of samples to apply in the margin
    durPulses = endPulses - startPulses;
    switch marginType
        case 'marginRatio'
            margins = floor(durPulses * marginValue);
        case 'marginTime'
            margins = floor(marginValue * ySampleRate)*ones(numDetPulses, 1);
    end
    startPulsesInd = startPulses + margins;
    endPulsesInd = endPulses - 1 - margins; % -1 because the ending indices are not included in the pulse (open interval)

    % Find the coefficients. 
    % They are the mean value of the pulse, once the margins have been applied
    coeff = zeros(numDetPulses, 1);
    for k = 1:numDetPulses
        ind = startPulsesInd(k):endPulsesInd(k);
        coeff(k) = mean(y(ind));
    end

end

function y = rectangularFilter(x, rectangleLimits, lags)

lengthX = size(x, 1);
lengthMask = max(rectangleLimits) - min(rectangleLimits) + 1;

if nargin < 3
minlag = -lengthMask - 1;
maxlag = lengthX - 1;
lags = minlag:maxlag;
end

acum = cumsum(x);
clear('x')
acum = [0; acum];

y = zeros(numel(lags), 1);
for lag_ind = 1:numel(lags)
    newLimits = rectangleLimits + lags(lag_ind);
    newLimits(newLimits < 1) = 0;
    newLimits(newLimits > lengthX) = lengthX;
    newLimits = newLimits + 1;
    
    y(lag_ind) = sum(acum(newLimits(:, 2))) - sum(acum(newLimits(:, 1)));
end

end

function shift = findShift(x, maskLimits, firstStep, minLag, maxLag)
% maskLimits. (numRectangles x 2) matrix. Each row are the first and last
% sample of a rectangle.

step = firstStep;

highLimit = maxLag;
lowLimit = minLag;

firstIteration = true;
while step >= 1
%     lags = maxLag:-step:minLag;
    lags = highLimit:-step:lowLimit;

    corr = rectangularFilter(x, maskLimits, lags);    
    corr = corr/max(corr); % Normalize the correlation values
    
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
    if isempty(ind)
        shift = 0;
        return;
    elseif ind == 1 || ind == length(lags)
        break;
    end
    
    lowLimit = lags(ind + 1);
    highLimit = lags(ind - 1);
    
    step = floor((highLimit - lowLimit)*0.1);
end

shift = lags(ind);

end
