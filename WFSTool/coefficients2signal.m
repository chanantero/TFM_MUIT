function [x, startInd, endInd] = coefficients2signal( coefficients, frequency, SampleRate )
% coefficients. (numChannels x numFreq)
% frequency. numFreq-element vector

numChann = size(coefficients, 1);
numFreq = size(coefficients, 2);

% First, reproduce in each channel and silence the rest
amplitude = abs(coefficients);
phase = angle(coefficients);
soundSamples = ceil(1*SampleRate);
silenceSamples = ceil(1*SampleRate);
[x_pre, startInd_pre, endInd_pre] = successiveChannelSinusoids( amplitude, phase, frequency, SampleRate, soundSamples, silenceSamples );
numPre = size(x_pre, 1);

% Then, reproduce everything at the same time
numSamples_main = ceil(SampleRate * 1); % 1 seconds
t0 = numPre/SampleRate;
t = t0 + (0:numSamples_main-1)'/SampleRate;
x_main = zeros(numSamples_main, numChann, numFreq);
for c = 1:numChann
    for f = 1:numFreq
        x_main(:, c, f) = amplitude(c, f) * cos(2*pi*frequency(f)*t + phase(c, f));
    end
end
x_main = sum(x_main, 3);

x = [x_pre; x_main];

% Set the time vector
startInd = [startInd_pre; numPre + 1];
endInd = [endInd_pre; numPre + numSamples_main];
end