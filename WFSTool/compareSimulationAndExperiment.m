function rat = compareSimulationAndExperiment( freqs, xOriginCoef, x, xPulseLimits, xSampleRate, y, ySampleRate, simul )
% Compare the simulation results with the experimental results.
% x. Reproduced signal. (numSamplesX x numChannelsX) matrix.
% xOriginCoef. Complex coefficients from which x was synthesized. (numPulses x numChannelsX x
% numFrequencies) matrix.
% freqs. Frequency of each sinusoidal signal. numFrequencies-element vector
% y. Received signal. (numSamplesY x numChannelsY) matrix
% simul. (numChannelsX x numFrequencies x numChannelsY).
% x is a signal with (numFrequencies + 1) parts. The first numFrequencies
% ones consist each one in a succession of sinusoidal signals for each
% channel while the rest of channels is silent. The last one is the
% complete signal on all channels.


% Get a pulse coefficient matrix of y that corresponds to the coefficient
% matrix of x
yPulseCoefMat = signal2pulseCoefficientMatrix(freq, xPulseCoefMat, xPulseLimits, xSampleRate, y, ySampleRate);

%  Identify the coefficients detected for each y channel, each x channel
%  (loudspeaker) and each frequency.
yCoef = zeros(numChannelsX, numChannelsY, numFrequencies);
for cy = 1:numChannelsY
    % Pulses corresponding to the f-th frequency
    
    % Each of those pulses correspond to one channel of x (loudspeaker)
    
end

% Chech if the total signal corresponds to the sum of the rest

% Relation of the received coefficients and the simulated ones.
rat = y_coef./simul;


end