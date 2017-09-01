function [expAcPath, freq] = getAcousticPath( xPulseLimits, freq, xPulseCoefMat, y, ySampleRate)
% Get the experimental acoustic path.
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

numPulses = size(xPulseCoefMat, 1);
numChannelsX = size(xPulseCoefMat, 2);
numChannelsY = size(y, 2);

% If there are repeated frequencies, combine them
[freq, ~, ic] = unique(freq);
numFrequencies = numel(freq);
xPulseCoefMat_uniqueFreq = zeros(numPulses, numChannelsX, numFrequencies);
for f = 1:numFrequencies
    xPulseCoefMat_uniqueFreq(:, :, f) = sum(xPulseCoefMat(:, :, ic == f), 3);
end

% Get a pulse coefficient matrix of y that corresponds to the coefficient
% matrix of x
yPulseCoefMat = signal2pulseCoefficientMatrix(xPulseLimits, freq, xPulseCoefMat_uniqueFreq, y, ySampleRate);

% Solve a linear system of equations. This way doesn't work if there are two
% frequencies that are the same and there is some pulse where those two
% "frequencies" (which are the same) are reproducing at the same time.
expAcPath = zeros(numChannelsY, numChannelsX, numFrequencies);
for f = 1:numFrequencies
    expAcPath(:, :, f) = (xPulseCoefMat_uniqueFreq\yPulseCoefMat).';
end

end

%% Old version
% function expAcPath = getAcousticPath( xPulseLimits, freq, xPulseCoefMat, y, ySampleRate)
% % Get the experimental acoustic path.
% % freq. numFrequencies-element vector. The i-th element is the i-th
% % frequency in Hz.
% % xPulseCoefMat. (numPulses x numChannelsX x numFrequencies) array. The
% % (i, j, k)-th element contains the complex coefficient of the i-th pulse
% % transmitted through the j-th channel in the k-th frequency.
% % xPulseLimits. (numPulses x 2) matrix. The i-th row contains the starting
% % and ending sample indices in which the i-th pulse occurs.
% % xSampleRate. Scalar. Sample rate of the signal x.
% % pulseInd. numSingularPulses-element vector. It might be that some of the pulses in x only
% % were transmitted through one channel and one frequency, in order to
% % perform some acoustic path measure. Those pulses are then associated to
% % one channel and one frequency. The indexes of those pulses are contained
% % in pulseInd. It has as many elements as pulses are singular (one
% % frequency and channel), i.e., intended for acoustic path measure. The
% % i-th element contains the index of the i-th singular pulse.
% % xChanInd. numSingularPulses-element vector. The i-th element contains the
% % channel associated with the i-th singular pulse.
% % freqInd. numSingularPulses-element vector. The i-th element contains the
% % frequency associated with the i-th singular pulse. 
% % y. (numSamplesY, numChannelsY) matrix. Received signal.
% % ySampleRate. Scalar. Sample rate of the received signal y.
% % simulCoef. (numChannelsX x numChannelsY x numFrequencies) array. The (i, j,
% % k)-th element contains the complex coefficient of the contribution that
% % procedes from the i-th x channel to the j-th y channel in the k-th
% % frequency.
% 
% numChannelsX = size(xPulseCoefMat, 2);
% numChannelsY = size(y, 2);
% numFrequencies = numel(freq);
% 
% % Get a pulse coefficient matrix of y that corresponds to the coefficient
% % matrix of x
% yPulseCoefMat = signal2pulseCoefficientMatrix(xPulseLimits, freq, xPulseCoefMat, y, ySampleRate);
% 
% %  Identify the coefficients detected for each y channel, each x channel
% %  (loudspeaker) and each frequency.
% [pulseInd, xChanInd, freqInd] = pulseCoefMat2singPulseInd( xPulseCoefMat );
% 
% yCoef = cell(numChannelsY, numChannelsX, numFrequencies);
% xCoef = cell(numChannelsX, numFrequencies);
% for p = 1:numel(pulseInd)
%     for cy = 1:numChannelsY
%         yCoef{cy, xChanInd(p), freqInd(p)} = [yCoef{cy, xChanInd(p), freqInd(p)}; yPulseCoefMat(pulseInd(p), cy, freqInd(p))];
%     end
%     xCoef{xChanInd(p), freqInd(p)} = [xCoef{xChanInd(p), freqInd(p)}; xPulseCoefMat(pulseInd(p), xChanInd(p), freqInd(p))];
% end
% 
% % Acoustic paths
% expAcPath = zeros(numChannelsY, numChannelsX, numFrequencies);
% for cx = 1:numChannelsX
%     for f = 1:numFrequencies
%         trans = xCoef{cx, f};
%         for cy = 1:numChannelsY
%             expAcPath(cy, cx, f) = mean(yCoef{cy, cx, f}./trans);
%         end
%     end
% end
% 
% % % Chech if the total signal corresponds to the sum of the rest
% % totalSignalCoef = yPulseCoefMat(end, :, :); % Last pulse
% % totalTheoricSignalCoef = sum(repmat(permute(xPulseCoefMat(end, :, :), [2, 1, 3]), [1, numChannelsY, 1]).*expAcPath, 1); % Sum of all channels and frequencies
% % totalSignalCoef./totalTheoricSignalCoef
% 
% end