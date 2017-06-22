function [ values, pulseLimits ] = pulseSignalParameters( x )
% x is a signal with complex values in form of pulses
N = size(x, 1);

% Calculate the median (mediana)
xSort = sort(abs(x));
median = xSort(floor(N/2 + 1));
threshold = median;

% Apply
active = abs(x) > threshold;
indActive = find(active);
diffActive = diff(indActive);
sep = find(diffActive > 1);

startPulse = [indActive(1); indActive(sep + 1)];
endPulse = [indActive(sep); indActive(end)];
numPulses = numel(startPulse);

values = zeros(numPulses, 1);
for k = 1:numPulses
    values(k) = mean(x(startPulse(k):endPulse(k)));
end

pulseLimits = [startPulse, endPulse];

end

