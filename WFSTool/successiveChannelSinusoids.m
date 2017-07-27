% % Vectorized version
% function [pulseCoefMat, pulseLimits] = successiveChannelSinusoids( coefMat, frequency, soundSamples, silenceSamples, numRep )
% % amplitude. (numChannels x numFreqs)
% 
% numChann = size(coefMat, 1);
% numFreq = numel(frequency);
% 
% % First, reproduce in each channel and silence the rest
% samplesPerChannel = soundSamples + silenceSamples;
% numPulses = numChann * numFreq * numRep;
% 
% pulseCoefMat_noSilences = zeros(numPulses, numChann, numFreq);
% pulseLimits = zeros(numPulses, 1);
% % pulseDurations = zeros(numPulses, 1);
% 
% [Ch, Fr, ~] = ndgrid(1:numChann, 1:numFreq, 1:numRep);
% c = Ch(:);
% f = Fr(:);
% for p = 1:numPulses     
%     startPulseSample = (p - 1)*samplesPerChannel + 1;
%     endPulseSample = (p - 1)*samplesPerChannel + soundSamples + 1;
%     
%     pulseLimits(p, 1) = startPulseSample;
%     pulseLimits(p, 2) = endPulseSample;
%     
% %     pulseDurations(p) = soundSamples;
%     
%     pulseCoefMat_noSilences(p, c(p), f(p)) = coefMat(c(p), f(p));
% end
% 
% % pulseDurations = reshape([soundSamples*ones(1, numPulses); silenceSamples*ones(1, numPulses)], [numPulses*2, 1]);
% 
% pulseLimits = [reshape(pulseLimits', [2*numPulses, 1]); pulseLimits(end) + silenceSamples];
% pulseCoefMat = zeros(numPulses*2, numChann, numFreq);
% pulseCoefMat(1:2:end, :, :) = pulseCoefMat_noSilences;
% 
% end

% Other version
function [pulseCoefMat, pulseLimits] = successiveChannelSinusoids( coefMat, frequency, soundSamples, silenceSamples, numRep )
% amplitude. (numChannels x numFreqs)

numChann = size(coefMat, 1);
numFreq = numel(frequency);

% First, reproduce in each channel and silence the rest
samplesPerChannel = soundSamples + silenceSamples;
numPulses = numChann * numFreq * numRep;

pulseCoefMat = zeros(numPulses, numChann, numFreq);

[Ch, Fr, ~] = ndgrid(1:numChann, 1:numFreq, 1:numRep);
c = Ch(:);
f = Fr(:);
for p = 1:numPulses        
    pulseCoefMat(p, c(p), f(p)) = coefMat(c(p), f(p));
end

% signal = pulseCoefMat2signal(frequency, coefMat, pulseLimits, SampleRate, startSample, endSample, false);

% % Set the markers vector
startPulseInd = (0:numPulses-1)'*samplesPerChannel;
endPulseInd = (0:numPulses-1)'*samplesPerChannel + soundSamples;

pulseLimits = [startPulseInd, endPulseInd];
end


% % Old version
% function [signal, startPulseInd, endPulseInd] = successiveChannelSinusoids_inds( amplitude, phase, frequency, SampleRate, soundSamples, silenceSamples, startSample, endSample )
% 
% numChann = size(amplitude, 1);
% numFreq = numel(frequency);
% 
% % First, reproduce in each channel and silence the rest
% amp = amplitude(:);
% phase = phase(:);
% freq = kron(frequency(:), ones(numChann, 1));
% 
% samplesPerChannel = soundSamples + silenceSamples;
% numCicles = numChann * numFreq;
% totalSamples = numCicles * samplesPerChannel;
% t = (0:totalSamples - 1)'/SampleRate;
% 
% endSample = min(endSample, totalSamples);
% startSample = max(startSample, 1);
% 
% indSamples = startSample:endSample;
% cicles = ceil(indSamples/samplesPerChannel);
% sampleInCicle = rem(indSamples, samplesPerChannel);
% sampleInCicle(sampleInCicle == 0) = samplesPerChannel;
% 
% signal = zeros(numel(indSamples), numChann);
% cicleIndices = unique(cicles);
% for k = cicleIndices
%     
%     c = rem(k, numChann);
%     if c == 0; c = numChann; end
%     
%     absInd = (k - 1)*samplesPerChannel + (1:soundSamples)';
%     fullCicle = zeros(samplesPerChannel, numChann);
%     fullCicle(1:soundSamples, c) = amp(k) * cos(2*pi*freq(k)*t(absInd) + phase(k)); 
%     
%     ind = find(cicles == k);
%     signal(ind, :) = fullCicle(sampleInCicle(ind), :);
%     
% end
% 
% % Set the markers vector. Starting indices of pulses
% startPulseIndTotal = (0:numCicles-1)'*samplesPerChannel + 1;
% endPulseIndTotal = (0:numCicles-1)'*samplesPerChannel + soundSamples;
% 
% startPulseInd = startPulseIndTotal(ismember(startPulseIndTotal, indSamples));% - startSample(1) + 1;
% endPulseInd = endPulseIndTotal(ismember(endPulseIndTotal, indSamples));% - startSample(1) + 1;
% 
% end