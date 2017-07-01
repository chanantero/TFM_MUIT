function [signal, startPulseInd, endPulseInd] = successiveChannelSinusoids_inds( amplitude, phase, frequency, SampleRate, soundSamples, silenceSamples, startSample, endSample )

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

endSample = min(endSample, totalSamples);
startSample = max(startSample, 1);

indSamples = startSample:endSample;
cicles = ceil(indSamples/samplesPerChannel);
sampleInCicle = rem(indSamples, samplesPerChannel);
sampleInCicle(sampleInCicle == 0) = samplesPerChannel;

signal = zeros(numel(indSamples), numChann);
cicleIndices = unique(cicles);
for k = cicleIndices
    
    c = rem(k, numChann);
    if c == 0; c = numChann; end
    
    absInd = (k - 1)*samplesPerChannel + (1:soundSamples)';
    fullCicle = zeros(samplesPerChannel, numChann);
    fullCicle(1:soundSamples, c) = amp(k) * cos(2*pi*freq(k)*t(absInd) + phase(k)); 
    
    ind = find(cicles == k);
    signal(ind, :) = fullCicle(sampleInCicle(ind), :);
    
end

% Set the markers vector. Starting indices of pulses
startPulseIndTotal = (0:numCicles-1)'*samplesPerChannel + 1;
endPulseIndTotal = (0:numCicles-1)'*samplesPerChannel + soundSamples;

startPulseInd = startPulseIndTotal(ismember(startPulseIndTotal, indSamples));% - startSample(1) + 1;
endPulseInd = endPulseIndTotal(ismember(endPulseIndTotal, indSamples));% - startSample(1) + 1;

end