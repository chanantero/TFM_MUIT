function [ output ] = applyDelay( input, Fs, delays )
% Input arguments:
% - input. Input signal. NxC array where N is the number of samples and C
% is the number of channels.
% - Fs. Scalar. Sampling frequency of 'input'.
% - delays. Delays. NxC array. The (i,j)-th element is the delay that has to
% be applied to the j-th channel at the time the i-th sample was supposed
% to be played if there was no delay. In other words, a delay of 2 seconds
% for the i-th sample means that at that time we must reproduce what should
% have been reproduced 2 seconds before. In other words, it is a delay from
% the point of view of the receiver, not the transmitter: received(t) =
% transmitted(t - delay(t))

[numSamples, numChannels] = size(input);

% Delays in terms of samples
delaySamples = round(delays * Fs);

% Delayed indices
indices = repmat((1:numSamples)', 1, numChannels); % Reference indices
indices = indices - delaySamples; % Delayed indices

% Apply delays
valid = indices > 0;
output = zeros(size(input));
output(valid) = input(indices(valid)); % Indices previous to the beginning of the signal are not set (they conserve the value 0)

end

