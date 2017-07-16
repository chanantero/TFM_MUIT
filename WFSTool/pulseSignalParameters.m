function [ values, pulseLimits ] = pulseSignalParameters( x, filterWidth, minDur, minSep, marginRatio )
% x is a signal with complex values in form of pulses

% Parse inputs
if nargin < 5
    marginRatio = 0.1;
end
if nargin < 4
    minSep = 1000;
end
if nargin < 3
    minDur = 1000;
end
if nargin < 2
    filterWidth = 600;
end

% Filter the signal. Low pass filter.
xFilt = filter(ones(filterWidth, 1), filterWidth, x); 

% Sort the values
xSort = sort(abs(xFilt));
% Gradient
grad = diff(xSort);
% Normalize the gradient
gradNorm = grad*numel(xSort)/max(xSort);
% Filter
gradFiltWidth = round(0.01*numel(gradNorm));
gradNormFilt = filter(ones(gradFiltWidth, 1), gradFiltWidth, gradNorm);

% Find local maxima. Each of them are divisions between levels of peaks
[maxtab, ~] = peakdet(gradNormFilt, 10);

% Find the mean value of first level
firstPeak = maxtab(1,1);
if size(maxtab, 1) < 2
    secondPeak = numel(xSort);
else
    secondPeak = maxtab(2, 1);
end
firstStepLevel = mean(xSort(firstPeak:secondPeak));
threshold = firstStepLevel/2;

% Apply threshold
active = abs(xFilt) > threshold;
indActive = find(active);
diffActive = diff(indActive);
sep = find(diffActive > 1);

% Join the pulses that are closer than the minimum separation
underSep = sep < minSep;
sep(underSep) = [];

startPulse = [indActive(1); indActive(sep + 1)];
endPulse = [indActive(sep); indActive(end)];

% Delete the pulses that are shorter than the minimum duration
dur = endPulse - startPulse;
underDur = dur < minDur;
startPulse(underDur) = [];
endPulse(underDur) = [];

% Get the complex coefficient of the pulse
numPulses = numel(startPulse);
values = zeros(numPulses, 1);
for k = 1:numPulses
    % Discard the extremes of the pulse
    durPulse = endPulse(k) - startPulse(k) + 1;
    margin = floor(durPulse * marginRatio);
    ind = (startPulse(k) + margin):(endPulse(k) - margin);
    values(k) = mean(x(ind));
end

pulseLimits = [startPulse, endPulse];

end


% % Calculate the median (mediana)
% % N = size(x, 1);
% % xSort = sort(abs(x));
% % median = xSort(floor(N/2 + 1));
% % threshold = median;
% threshold = mean(abs(x));
% 
% % Apply
% active = abs(x) > threshold;
% indActive = find(active);
% diffActive = diff(indActive);
% sep = find(diffActive > 1);
% 
% startPulse = [indActive(1); indActive(sep + 1)];
% endPulse = [indActive(sep); indActive(end)];
% numPulses = numel(startPulse);
% 
% values = zeros(numPulses, 1);
% for k = 1:numPulses
%     values(k) = mean(x(startPulse(k):endPulse(k)));
% end
% 
% pulseLimits = [startPulse, endPulse];
