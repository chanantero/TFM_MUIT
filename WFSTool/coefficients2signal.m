function [x, outOfRange, pulseCoefMat, pulseLimits ] = coefficients2signal( coefficients, frequency, SampleRate, startSample, endSample, onlyPulseInfoFlag )
% coefficients. (numChannels x numFreq)
% frequency. numFreq-element vector
if nargin < 6
    onlyPulseInfoFlag = false;
end


[ pulseCoefMat, pulseLimits ] = coefficients2pulseSignalParameters( coefficients, frequency, SampleRate, 'prelude' );

if nargin < 4
    startSample = pulseLimits(1, 1);
    endSample = pulseLimits(end, 2);
end

if ~onlyPulseInfoFlag
    [x, outOfRange] = pulseCoefMat2signal(frequency, pulseCoefMat, pulseLimits, SampleRate, startSample, endSample, 'sample'); % endSample+1 because of the way pulseCoefMat2signal works
else
    x = [];
    outOfRange = [];
end

%% Singular Pulse information
% % sPMat should be the output in case we wanted to use that information.
% [Ch, Fr, ~] = ndgrid(1:numChann, 1:numFreq, 1:numRep);
% numPrePulses = numel(Ch);
% sPind = (1:numPrePulses)'; % Singular pulse indices
% sPChInd = Ch(:); % Singular pulse channel indices
% sPFreqInd = Fr(:); % Singular pulse frequency indices
% sPMat = [sPind, sPChInd, sPFreqInd]; % Singular pulse information matrix

end