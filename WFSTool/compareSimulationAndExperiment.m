function rat = compareSimulationAndExperiment( freqs, xOriginCoef, x, x_limitSamples, xSampleRate, y, ySampleRate, simul )
% Compare the simulation results with the experimental results.
% x. Reproduced signal. (numSamplesX x numChannelsX) matrix.
% xOriginCoef. Complex coefficients from which x was synthesized. (numChannelsX x
% numFrequencies) matrix.
% freqs. Frequency of each sinusoidal signal. numFrequencies-element vector
% y. Received signal. (numSamplesY x numChannelsY) matrix
% simul. (numChannelsX x numFrequencies x numChannelsY).
% x is a signal with (numFrequencies + 1) parts. The first numFrequencies
% ones consist each one in a succession of sinusoidal signals for each
% channel while the rest of channels is silent. The last one is the
% complete signal on all channels.

numChannelsX = size(x, 2); % Number of channels of the reproduction device
numChannelsY = size(y, 2); % Number of channels of the recorder device
numFrequencies = numel(freq);

% Pulses that should be detected in y. This is, pulses that are in x but
% counting the contribution of all x channels together
xGroupCoef = cell(numFrequencies, 1);
xGroupPulseLimits = cell(numFrequencies, 1);
for f_ind = 1:numFrequencies
    
end

tol = 0.1; % Tolerance
[corrInd, y_coef, y_pulseLimits] = associateSinPulses(freqs, xPulseLimits, xSampleRate, y, ySampleRate, tol);

% Map the corresponding indices into the original x signal pulse structure


divInd;

y_coef = zeros(numChannelsX, numFrequencies, numChannelsY);
for f = 1:numFrequencies
    y_part = y(divInd(f):divInd(f+1), :);
    iq = real2IQ(y_part, y_SampleRate, freq(f));
        
    for r = 1:numChannelsY
        [ values, pulseLimits ] = pulseSignalParameters( iq(:, r) );
        % We expect numChannels pulses, and hence, values
        % Divide each real value by it's source coefficient
        if numel(values) == numChannelsX
            y_coef(:, f, r) = values.';
        else
            error('compareSimulationAndExperiment:WrongSignal', 'Wrong number of pulses detected')
        end
    end
    
end

% Chech if the total signal corresponds to the sum of the rest

% Relation of the received coefficients and the simulated ones.
rat = y_coef./simul;


end