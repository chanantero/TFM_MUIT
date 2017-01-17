function [ backDelay ] = forward2BackwardDelay( forwDelay, Fs, firstSample, lastSample )
% Input arguments:
% - forwDelay. Forward delay. Vector of N elements, each one corresponding to
% one sample.
% - Fs. Sampling frequency. Scalar.
% - firstSample. First sample in which the backward delay will be
% interpolated. It is referenced from the first sample of forwDelay.
% - lastSample. Last sample in which the backward delay will be
% interpolated. It is referenced from the first sample of forwDelay.
% 
% Output arguments:
% - backDelay. Backward delay. Vector of N elements, each one corresponding
% to one sample.

forwardTime = (0:numel(forwDelay)-1)'/Fs;
backwardTime = forwardTime + forwDelay;

queryBackwardTime = (firstSample:lastSample)/Fs;

backDelay = interp1(backwardTime, forwDelay, queryBackwardTime, 'nearest');

end

