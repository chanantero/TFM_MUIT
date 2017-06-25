function yPulseCoefMat = signal2pulseCoefficientMatrix(freq, xPulseCoefMat, xPulseLimits, xSampleRate, y, ySampleRate)
% freq. Frequencies
% xPulseCoefMat. Original pulse coefficient matrix. (numPulsesX x
% numChannelsX x 
% y. Received signal.

numChannelsY = size(y, 2); % Number of channels of the recorder device
numFrequencies = numel(freq);

% Pulses that should be detected in y. This is, pulses that are in x but
% counting the contribution of all x channels together
activePulses = zeros(numPulsesX, numFrequencies);
xGroupPulseLimits = cell(numFrequencies, 1);
for f_ind = 1:numFrequencies
    xCoefMat1freq = xPulseCoefMat(:, :, f_ind);
    activePulses1freq = any(xCoefMat1freq ~= 0, 2);
    activePulses(:, f_ind) = activePulses1freq;
    xGroupPulseLimits{f_ind} = xPulseLimits(activePulses1freq);
end

tol = 0.1; % Tolerance
[corrInd, y_coef, ~] = associateSinPulses(freqs, xGroupPulseLimits, xSampleRate, y, ySampleRate, tol);

% Map the corresponding indices into the original x signal pulse structure.
% Remember that the pulse limits in y can be very different from the ones
% in x for each channel and frequency
yPulseCoefMat = zeros(numPulsesX, numChannelsY, numFrequencies);
for f_ind = 1:numFrequencies
    for cy = 1:numChannelsY
        yPulseCoefMat(activePulses(:, f_ind), cy, f_ind) = y_coef{f_ind, cy}(corrInd{f_ind, cy});
    end
end

end