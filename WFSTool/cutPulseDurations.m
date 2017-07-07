function [newPulseDurations, indPulses] = cutPulseDurations(pulseDurations, startMarker, endMarker)
% pulseDurations. (numPulses)-element vector.

pulseLimits = [0; cumsum(pulseDurations)];

[newPulseLimits, indPulses] = cutPulseLimits(pulseLimits, startMarker, endMarker);

newPulseDurations = diff(newPulseLimits);

end