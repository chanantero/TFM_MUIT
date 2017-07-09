function [newIntervalLimits, indIntervals] = intervalSelection(intervalLimits, startMarker, endMarker)
% intervalLimits. (numIntervals x 2) matrix. The (i, 1)-th element is the
% starting point of the i-th interval. The (i, 2)-th element is the ending
% point of the i-th interval. The i-th interval exist in the interval
% [intervalLimits(i, 1), intervalLimits(i, 2)[ (notice the open interval at the
% end interval, this is, the ending point is not included in the interval.
% startMarker. Starting point of the output sequence. It is included in the
% sequence.
% endMarker. Ending point of the output sequence sequence. It is
% interval-open, this is, it is not included.
% Output arguments:
% - newIntervalLimits.
% - indIntervals. Interval indices included in the selection

numIntervals = size(intervalLimits, 1);

% Transform the pulseLimits matrix to a column vector. So, the silence
% between pulses is now considered as pulses of zero amplitude
intervLim_vec = reshape(intervalLimits', [numIntervals*2, 1]);

% firstInterval in [1, numIntervals + 1]
firstInterval = find(intervLim_vec - startMarker > 0, 1, 'first') - 1;
if isempty(firstInterval)
    firstInterval = numIntervals + 1;
end
if firstInterval == 0
    firstInterval = 1;
end

% lastInterval in [0, numIntervals]
lastInterval = find(endMarker - intervLim_vec > 0, 1, 'last');
if isempty(lastInterval)
    lastInterval = 0;
end
if lastInterval > numIntervals
    lastInterval = numIntervals;
end

indIntervals = firstInterval:lastInterval;

% Select only odd intervals, because they are the real ones. The even
% intervals are artificial null intermediate intervals.
indIntervals = indIntervals(rem(indIntervals, 2) == 1);
indIntervals = (indIntervals + 1)/2;

% Determine the new interval limits. The interval limits of all intervals
% except the first and the last one are the original ones. However, the
% first and last intervals can be cutted by the start and end markers.
newIntervalLimits = intervalLimits(indIntervals, :);
newIntervalLimits(1, 1) = max(startMarker, newIntervalLimits(1, 1));
newIntervalLimits(end, 2) = min(endMarker, newIntervalLimits(end, 2));


end
