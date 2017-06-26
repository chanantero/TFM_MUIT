function [acPath, rat] = compareSimulationAndExperiment( freq, xPulseCoefMat, xPulseLimits, xSampleRate, pulseInd, xChanInd, freqInd, y, ySampleRate, simul )
% Compare the simulation results with the experimental results.
% freq. numFrequencies-element vector. The i-th element is the i-th
% frequency in Hz.
% xPulseCoefMat. (numPulses x numChannelsX x numFrequencies) array. The
% (i, j, k)-th element contains the complex coefficient of the i-th pulse
% transmitted through the j-th channel in the k-th frequency.
% xPulseLimits. (numPulses x 2) matrix. The i-th row contains the starting
% and ending sample indices in which the i-th pulse occurs.
% xSampleRate. Scalar. Sample rate of the signal x.
% pulseInd. numSingularPulses-element vector. It might be that some of the pulses in x only
% were transmitted through one channel and one frequency, in order to
% perform some acoustic path measure. Those pulses are then associated to
% one channel and one frequency. The indexes of those pulses are contained
% in pulseInd. It has as many elements as pulses are singular (one
% frequency and channel), i.e., intended for acoustic path measure. The
% i-th element contains the index of the i-th singular pulse.
% xChanInd. numSingularPulses-element vector. The i-th element contains the
% channel associated with the i-th singular pulse.
% freqInd. numSingularPulses-element vector. The i-th element contains the
% frequency associated with the i-th singular pulse. 
% y. (numSamplesY, numChannelsY) matrix. Received signal.
% ySampleRate. Scalar. Sample rate of the received signal y.
% simul. (numChannelsX x numChannelsY x numFrequencies) array. The (i, j,
% k)-th element contains the complex coefficient of the contribution that
% procedes from the i-th x channel to the j-th y channel in the k-th
% frequency.

% Get a pulse coefficient matrix of y that corresponds to the coefficient
% matrix of x
yPulseCoefMat = signal2pulseCoefficientMatrix(freq, xPulseCoefMat, xPulseLimits, xSampleRate, y, ySampleRate);

%  Identify the coefficients detected for each y channel, each x channel
%  (loudspeaker) and each frequency.
yCoef = zeros(numChannelsX, numChannelsY, numFrequencies);
xCoef = zeros(numChannelsX, numFrequencies);
for p = 1:numel(pulseInd)
    yCoef(xChanInd(p), :, freqInd(p)) = yPulseCoefMat(pulseInd(p), :, freqInd(p));
    xCoef(xChanInd(p), freqInd(p)) = xPulseCoefMat(pulseInd(p), xChanInd(p), freqInd(p));
end

% Acoustic paths
acPath = zeros(numChannelsX, numChannelsY, numFrequencies);
for cy = 1:numChanY
    acPath(:, cy, :) = yCoef(:, cy, :)./permute(xCoef, [1, 3, 2]);
end

% Chech if the total signal corresponds to the sum of the rest
totalSignalCoef = yPulseCoefMat(end, :, :); % Last pulse
totalTheoricSignalCoef = sum(repmat(permute(xPulseCoefMat(end, :, :), [2, 1, 3]), [1, numChannelsY, 1]).*acPath, 1); % Sum of all channels and frequencies

% Relation of the received coefficients and the simulated ones.
rat = yCoef./simul;

end