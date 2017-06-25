function [corrInd, y_coef, y_pulseLimits] = associateSinPulses(freqs, xPulseLimits, xSampleRate, y, ySampleRate, durationTolerance)

% y. Received signal. (numRecSamples x numChannels) matrix
% xPulseLimits. Sample limits of the pulses in the original signal x. 
% numFrequencies-element cell column vector. The i-th cell contains the sample
% limits for the i-th frequency.

numFrequencies = numel(freqs);
numChannelsY = size(y, 2);

% Detect the pulses in y
y_coef = cell(numFrequencies, numChannelsY);
y_pulseLimits = cell(numFrequencies, numChannelsY);
for f_ind = 1:numFrequencies
    f = freqs(f_ind);
    iq = real2IQ(y, ySampleRate, f, 5);
    
    for cy = 1:numChannelsY
        [ y_coef{f_ind, cy}, y_pulseLimits{f_ind, cy} ] = pulseSignalParameters( iq(:, cy) );
    end
end

% Association of pulses
corrInd = cell(numFrequencies, numChannelsY);
for f_ind = 1:numFrequencies
    xPulsLim = xPulseLimits{f_ind};
    durPulsesX = diff(xPulsLim, 2)/xSampleRate;
    marg = durPulsesX * durationTolerance;
        
    for cy = 1:numChannelsY
        yPulsLim = y_pulseLimits{f_ind, cy};
        durPulsesY = diff(yPulsLim, 2)/y_SampleRate;        
        corrInd{f_ind, cy} = associateValues(durPulsesX, durPulsesY, marg);        
    end
end

end


function corrInd = associateValues(x, y, margins)
% x. Column vector with original values
% y. Column vector with distorded values.
numValuesX = numel(x);
numValuesY = numel(y);

corrInd = zeros(numValuesX, 1);

yi = 1; % y index
for xi = 1:numValuesX
    if abs(y(yi) - x(xi)) <= margins(xi)
        corrInd(xi) = yi;
        yi = yi + 1;
        if yi > numValuesY
            break;
        end
    end
end

end
