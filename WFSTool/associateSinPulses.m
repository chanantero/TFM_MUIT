function [corrInd, y_coef, y_pulseSampleLimits] = associateSinPulses(freqs, xPulseLimits, y, ySampleRate)

% y. Received signal. (numRecSamples x numChannels) matrix
% xPulseLimits. Time limits of the pulses in the original signal x. 
% numFrequencies-element cell column vector. The i-th cell contains the sample
% limits for the i-th frequency.

numFrequencies = numel(freqs);
numChannelsY = size(y, 2);

% Detect the pulses in y
y_coef = cell(numFrequencies, numChannelsY);
corrInd = cell(numFrequencies, numChannelsY);
y_pulseSampleLimits = cell(numFrequencies, numChannelsY);
for f_ind = 1:numFrequencies
    f = freqs(f_ind);
    iq = real2IQ(y, ySampleRate, f, 15);
    
    for cy = 1:numChannelsY
        [y_coef{f_ind, cy}, corrInd{f_ind, cy}, y_pulseSampleLimits{f_ind, cy}] = detectPulseSignal(xPulseLimits{f_ind}, iq(:, cy), ySampleRate, 'marginTime', 0.1);
    end
     
%     for cy = 1:numChannelsY
%         [ y_coef{f_ind, cy}, y_pulseSampleLimits{f_ind, cy} ] = pulseSignalParameters( iq(:, cy) );
%     end
end

end


%% Old version.
% It doesn't use the function detectPulseSignal, but uses
% pulseSignalParameters and associateValues
% function [corrInd, y_coef, y_pulseSampleLimits] = associateSinPulses(freqs, xPulseLimits, y, ySampleRate)
% 
% % y. Received signal. (numRecSamples x numChannels) matrix
% % xPulseLimits. Time limits of the pulses in the original signal x. 
% % numFrequencies-element cell column vector. The i-th cell contains the sample
% % limits for the i-th frequency.
% 
% numFrequencies = numel(freqs);
% numChannelsY = size(y, 2);
% 
% % Detect the pulses in y
% y_coef = cell(numFrequencies, numChannelsY);
% y_pulseSampleLimits = cell(numFrequencies, numChannelsY);
% for f_ind = 1:numFrequencies
%     f = freqs(f_ind);
%     iq = real2IQ(y, ySampleRate, f, 500);
%     
%     % For Debugging
%     % Recover the signal
% %     aux = WFSanalyzer;
% %     aux.analyzeSignal(y, ySampleRate);
% 
%     
%     for cy = 1:numChannelsY
%         [ y_coef{f_ind, cy}, y_pulseSampleLimits{f_ind, cy} ] = pulseSignalParameters( iq(:, cy) );
%     end
%     
%     [coeff, corrInd] = detectPulseSignal(xPulseLimits, iq, ySampleRate, 0.1);
% end
% 
% % Association of pulses
% corrInd = cell(numFrequencies, numChannelsY);
% for f_ind = 1:numFrequencies
%     xPulsLim = xPulseLimits{f_ind};
%     durPulsesX = diff(xPulsLim, 1, 2);
%         
%     for cy = 1:numChannelsY
%         yPulsSampleLim = y_pulseSampleLimits{f_ind, cy};
%         durPulsesY = diff(yPulsSampleLim, 1, 2)/ySampleRate;        
%         corrInd{f_ind, cy} = associateValues(durPulsesX, durPulsesY);        
%     end
% end
% 
% end
% 
% 
% function corrInd = associateValues(x, y)
% % x. Column vector with original values
% % y. Column vector with distorded values.
% 
% numValuesX = numel(x);
% 
% [ insertionCost, deletionCost, substitutionCost ] = genCostMatrices( x, y );
% alignment = NeedlemanWunsch(insertionCost, deletionCost, substitutionCost);
% 
% corrInd = zeros(numValuesX, 1);
% x_ind = 1;
% y_ind = 1;
% for k = 1:numel(alignment)
%     if alignment(k) == 'substitute'
%         corrInd(x_ind) = y_ind;
%         x_ind = x_ind + 1;
%         y_ind = y_ind + 1;
%     elseif alignment(k) == 'delete'
%         x_ind = x_ind + 1;
%     elseif alignment(k) == 'insert'
%         y_ind = y_ind + 1;
%     end   
% end
% 
% % % Old way
% % % margin is an input
% % numValuesX = numel(x);
% % numValuesY = numel(y);
% % 
% % corrInd = zeros(numValuesX, 1);
% % 
% % yi = 1; % y index
% % for xi = 1:numValuesX
% %     if abs(y(yi) - x(xi)) <= margins(xi)
% %         corrInd(xi) = yi;
% %         yi = yi + 1;
% %         if yi > numValuesY
% %             break;
% %         end
% %     end
% % end
% 
% end
