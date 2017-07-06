function [x, pulseCoefMat, pulseLimits, sPMat ] = coefficients2signal( coefficients, frequency, SampleRate, startSample, endSample, onlyPulseInfoFlag )
% coefficients. (numChannels x numFreq)
% frequency. numFreq-element vector
if nargin < 6
    onlyPulseInfoFlag = false;
end
if nargin < 4
    startSample = 1;
    endSample = Inf;
end

numChann = size(coefficients, 1);
numFreq = size(coefficients, 2);

% Prelude parameters
soundTime = 1;
silenceTime = 1;
numRep = 2;

% Main signal parameters
mainTime = 2;

soundSamples = ceil(soundTime * SampleRate);
silenceSamples = ceil(silenceTime * SampleRate);
mainSamples = ceil(mainTime*SampleRate);

% First, reproduce in each channel and silence the rest
[pulseCoefMat_pre, pulseLimitsPre] = successiveChannelSinusoids( coefficients, frequency, soundSamples, silenceSamples, numRep);

% Then, reproduce everything at the same time
pulseCoefMat_main = permute(coefficients, [3, 1, 2]);
pulseLimits = [pulseLimitsPre; pulseLimitsPre(end) + mainSamples];
pulseCoefMat = [pulseCoefMat_pre; pulseCoefMat_main];

if ~onlyPulseInfoFlag
    x = pulseCoefMat2signal(frequency, pulseCoefMat, pulseLimits, SampleRate, startSample, endSample+1, 'sample'); % endSample+1 because of the way pulseCoefMat2signal works
else
    x = [];
end

%% Pulse information

% Create the pulse coefficient matrix
[Ch, Fr, ~] = ndgrid(1:numChann, 1:numFreq, 1:numRep);
numPrePulses = numel(Ch);
sPind = (1:numPrePulses)'; % Singular pulse indices
sPChInd = Ch(:); % Singular pulse channel indices
sPFreqInd = Fr(:); % Singular pulse frequency indices
sPMat = [sPind, sPChInd, sPFreqInd]; % Singular pulse information matrix

end