function [newPulseLimits, indPulses] = cutPulseLimits(pulseLimits, startMarker, endMarker)

% startMarker = max(0, startMarker);
% endMarker = min(pulseLimits(end), endMarker);
startLimits = pulseLimits(:, 1);
endLimits = pulseLimits(:, 2);

firstPulse = find(startMarker < endLimits, 1, 'first');
lastPulse = find(endMarker > startLimits, 1, 'last');
indPulses = firstPulse:lastPulse;

if ~isempty(indPulses)
    firstPulseStart = max(startMarker, startLimits(firstPulse));
    lastPulseEnd = min(endMarker, endLimits(lastPulse));
    
    newStartLimits = startLimits(indPulses);
    newStartLimits(1) = firstPulseStart;
    newEndLimits = endLimits(indPulses);
    newEndLimits(end) = lastPulseEnd;
    
    newPulseLimits = [newStartLimits, newEndLimits];
else
    newPulseLimits = zeros(0, 2);
end

end