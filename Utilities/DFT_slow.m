function result = DFT_slow(signals, sampleRate, frequencies)
% Input arguments:
% - signals. (numSamples x numSignals)
% - frequencies. (numFrequencies x 1)
% Output arguments:
% - result. (numFrequencies x numSignals)

% Transform the impulse response of the loudspeakers to frequency domain
[numSamples, numSignals] = size(signals);
numFrequencies = numel(frequencies);

t = (0:numSamples-1)'/sampleRate;

result = zeros(numFrequencies, numSignals);
for indF = 1:numFrequencies
    f = frequencies(indF);
    tone = exp(-1i*2*pi*f*t);
    result(indF, :) = sum(signals.*repmat(tone, [1, numSignals]), 1)/sampleRate;
end

end