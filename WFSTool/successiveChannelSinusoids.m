function [signal, startInd, endInd] = successiveChannelSinusoids( amplitude, phase, frequency, SampleRate, soundSamples, silenceSamples )

numChann = size(amplitude, 1);
numFreq = numel(frequency);

% First, reproduce in each channel and silence the rest
amp = amplitude(:);
phase = phase(:);
freq = kron(frequency(:), ones(numChann, 1));

samplesPerChannel = soundSamples + silenceSamples;
numCicles = numChann * numFreq;
totalSamples = numCicles * samplesPerChannel;
t = (0:totalSamples - 1)'/SampleRate;

signal = zeros(totalSamples, numChann);
for k = 1:numCicles
    c = rem(k, numChann);
    if c == 0; c = numChann; end
    ind = (k - 1)*samplesPerChannel + (1:soundSamples)';
    signal(ind, c) = amp(k) * cos(2*pi*freq(k)*t(ind) + phase(k));
end

% Set the markers vector. Starting indices of 
startInd = (0:numCicles-1)'*samplesPerChannel + 1;
endInd = (0:numCicles-1)'*samplesPerChannel + soundSamples;

end

% function [signal] = successiveChannelSinusoids( amplitude, phase, frequency, SampleRate, soundCicles, silenceCicles )
% % amplitude: numChannels-element vector
% % phase: numChannels-element vector
% % frequency. numFrequencies-element vector.
% % SampleRate. Scalar.
% % soundCicles. Scalar.
% % silenceCicles. Scalar.
% 
% numChannels = numel(amplitude);
% 
% cicleSamples = ceil(SampleRate/frequency*soundCicles);
% silenceSamples = ceil(SampleRate/frequency*silenceCicles);
% samplesPerChannel = cicleSamples + silenceSamples;
% totalSamples = numChannels*samplesPerChannel;
% 
% signal = zeros(totalSamples, numChannels);
% t = (0:totalSamples - 1)'/SampleRate;
% for k = 1:numChannels
%     ind = (k - 1)*samplesPerChannel + (1:cicleSamples);
%     signal(ind, k) = amplitude(k) * cos(2*pi*frequency*t(ind) + phase(k));
% end
% 
% % Determine the complex coefficients of the source signal for each channel
% % t_ini = samplesPerChannel*(0:numChannels - 1)/SampleRate;
% % phas = phase - 2*pi*frequency*t_ini;
% % complexCoef = amplitude*exp(1i*phas);
% 
% end