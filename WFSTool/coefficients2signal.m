function [x, outOfRange, pulseCoefMat, pulseLimits ] = coefficients2signal( coefficients, frequency, SampleRate, startSample, endSample, onlyPulseInfoFlag )
% coefficients. (numChannels x numFreq)
% frequency. numFreq-element vector
if nargin < 6
    onlyPulseInfoFlag = false;
end

numChann = size(coefficients, 1);
numFreq = size(coefficients, 2);

% Prelude parameters
soundTime = 1;
silenceTime = 1;
numRep = 5;

% Main signal parameters
mainTime = 5;

soundSamples = ceil(soundTime * SampleRate);
silenceSamples = ceil(silenceTime * SampleRate);
mainSamples = ceil(mainTime*SampleRate);

% First, reproduce in each channel and silence the rest
[pulseCoefMat_pre, pulseLimitsPre] = successiveChannelSinusoids( coefficients, frequency, soundSamples, silenceSamples, numRep);

% % Change amplitude of second repetition
% pulseCoefMat_pre(numChann*numFreq+1:end, :) = pulseCoefMat_pre(numChann*numFreq+1:end, :)*1.3;

% Then, reproduce everything at the same time
pulseCoefMat_main = permute(coefficients, [3, 1, 2]);
pulseLimits = [pulseLimitsPre; [pulseLimitsPre(end) + silenceSamples, pulseLimitsPre(end) + silenceSamples + mainSamples]; ...
    [pulseLimitsPre(end) + silenceSamples + mainSamples, pulseLimitsPre(end) + silenceSamples + mainSamples + silenceSamples]];
pulseCoefMat = [pulseCoefMat_pre; pulseCoefMat_main; zeros(1, numChann, numFreq)];


if nargin < 4
    startSample = pulseLimits(1, 1);
    endSample = pulseLimits(end, 2);
end

if ~onlyPulseInfoFlag
    [x, outOfRange] = pulseCoefMat2signal(frequency, pulseCoefMat, pulseLimits, SampleRate, startSample, endSample, 'sample'); % endSample+1 because of the way pulseCoefMat2signal works
else
    x = [];
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