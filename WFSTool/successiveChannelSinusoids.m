function [signal, complexCoef] = successiveChannelSinusoids( amplitude, phase, frequency, SampleRate, soundCicles, silenceCicles )
% amplitude: numChannels-element vector
% phase: numChannels-element vector
% frequency. Scalar.
% SampleRate. Scalar.
% soundCicles. Scalar.
% silenceCicles. Scalar.

numChannels = numel(amplitude);

cicleSamples = ceil(SampleRate/frequency*soundCicles);
silenceSamples = ceil(SampleRate/frequency*silenceCicles);
samplesPerChannel = cicleSamples + silenceSamples;

signal = zeros(numChannels*samplesPerChannel, numChannels);
for k = 1:numChannels
    ind = (k - 1)*samplesPerChannel + (1:cicleSamples);
    t = (0:cicleSamples - 1)/SampleRate;
    signal(ind, k) = amplitude(k) * cos(2*pi*frequency*t + phase(k));
end

% Determine the complex coefficients of the source signal for each channel
t_ini = samplesPerChannel*(0:numChannels - 1)/SampleRate;
phas = phase - 2*pi*frequency*t_ini;
complexCoef = amplitude*exp(1i*phas);

end