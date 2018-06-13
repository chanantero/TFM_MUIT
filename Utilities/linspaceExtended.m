function [ vec ] = linspaceExtended( positions, numIntervals)

positions = positions(:);
numIntervals = numIntervals(:);

numSegments = numel(positions) - 1;

numPoints = sum(numIntervals) + 1;
acum = cumsum(numIntervals);
startIntervalInd = [1; acum(1:end-1) + 1];
endIntervalInd = acum;
widthInterval = positions(2:end) - positions(1:end-1);
widthSubInterval = widthInterval./numIntervals;

vec = zeros(numPoints, 1);
for k = 1:numSegments
    vec(startIntervalInd(k):endIntervalInd(k)) =...
        positions(k):widthSubInterval(k):positions(k + 1) - widthSubInterval(k);
end
vec(end) = positions(end);

end

