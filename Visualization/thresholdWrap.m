function [ thresholds, changeDirections, allAboveThres, allUnderThres ] = thresholdWrap( points, func, thresDim, thresValue, sampVec, tolerance )
% Input arguments:
% - points. Cell vector with N components, beint N the number of
% dimensions. Each cell contains a numerical or categorical vector. A
% combination of all defines a N-dimensional grid with all the points for
% which a threshold is going to be calculated.
% - thresDim. Dimension on which the threshold search is going to be
% performed.
% - thresValue. Value of the threshold.
% - func. Function R^N --> R.
% - limits. Limits of the sampling segment along the threDim-th dimension.
% - numSamp. Number of points of the sampling segment.
% - tolerance. Precission that puts a stop to the shoot method.
% Aditional notes.
% The usual is that the thresDim-th cell of points has only one scalar,
% since the information that contains is not going to be used

points{thresDim} = sampVec;

samples = animation.calculateValues(points, {func});
samples = samples{1};

[ thresholds, changeDirections, allAboveThres, allUnderThres ] = searchThreshold( points, samples, func, thresDim, thresValue, tolerance );
end

