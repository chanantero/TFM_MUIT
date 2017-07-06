function [newPulseLimits, indPulses] = cutPulseLimits(pulseLimits, startMarker, endMarker)
% pulseLimits. (numPulses+1)-element vector.

numPulses = numel(pulseLimits) - 1;

firstPulse = find(pulseLimits - startMarker > 0, 1, 'first') - 1;
if isempty(firstPulse)
    firstPulse = numPulses + 1;
end
lastPulse = find(endMarker - pulseLimits > 0, 1, 'last');
if isempty(lastPulse)
    lastPulse = 0;
end

if firstPulse == 0
    firstPulse = 1;
end

indPulses = firstPulse:lastPulse;

if ~isempty(indPulses)    
    newPulseLimits = pulseLimits(firstPulse:lastPulse+1);
    newPulseLimits(1) = startMarker;
    newPulseLimits(end) = endMarker;
else
    newPulseLimits = zeros(0, 0);
end

end