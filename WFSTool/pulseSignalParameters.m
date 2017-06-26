function [ values, pulseLimits ] = pulseSignalParameters( x )
% x is a signal with complex values in form of pulses

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

% Other way
grad = diff(abs(x));
threshold_grad;
limits = find(grad > threshold_grad);

% Mean of the space between the limits
startPulse = [1; limits(1:end)];
endPulse = [limits(1:end); numel(x)];
numPulses = numel(startPulse);

values = zeros(numPulses, 1);
for k = 1:numPulses
    values(k) = mean(x(startPulse(k):endPulse(k)));
end

% Decide which of them are pulses and which of them are silence



end

