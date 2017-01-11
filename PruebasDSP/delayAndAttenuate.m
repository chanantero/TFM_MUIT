function [ output ] = delayAndAttenuate( input, Fs, delays, attenuations)
% Takes one-channel signal and applies delays and attenuations for each of
% the output channels
%
% Input arguments:
% - input. input signal, Nx1 numeric vector.
% - Fs. Sampling frequency.
% - delays. 1xC vector, where C is the number of channels of the output
% signal. The i-th element is the delay of the i-th channel
% - attenuation. 1xC vector, where C is the number of channels of the output
% signal. The i-th element is the attenuation of the i-th channel

numOutputChannels = numel(delays);
numSamples = numel(input);

% Apply the delay
delays = repmat(delays, numSamples, 1);
delayed = applyDelay( input, Fs, delays );

% Apply the attenuation
attenuations = repmat(attenuations, numSamples, 1);
output = delayed.*attenuations;

end

