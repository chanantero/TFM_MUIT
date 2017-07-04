function [sP, sPCh, sPFreq] = pulseCoefMat2singPulseInd( pulseCoefMat )
% pulseCoefMat. (numPulses x numChannels x numFrequencies)

% Find the pulses that only have one coefficient different from 0.
active = abs(pulseCoefMat) > 0;
activeInd = find(active);

numActiveCoef = sum(sum(active, 3), 2);
sP = find(numActiveCoef == 1);
numSingPulses = numel(sP);

[pulseInd, chanInd, freqInd] = ind2sub(size(pulseCoefMat), activeInd);

sPCh = zeros(numSingPulses, 1);
sPFreq = zeros(numSingPulses, 1);
for p = 1:numSingPulses
    ind = pulseInd == sP(p);
    sPCh(p) = chanInd(ind);
    sPFreq(p) = freqInd(ind);
end

end