function [x, startInd, endInd] = coefficients2signal_inds( coefficients, frequency, SampleRate, startSample, endSample, onlyPulseInfoFlag )
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

amplitude = abs(coefficients);
phase = angle(coefficients);

% Prelude parameters
soundSamples = ceil(0.5*SampleRate);
silenceSamples = ceil(0.5*SampleRate);
samplesPerChannel = soundSamples + silenceSamples;
numCicles = numChann * numFreq;
totalPre = numCicles * samplesPerChannel;

% Main signal parameters
totalMain = ceil(1*SampleRate);

% Indices
startSamplePre = startSample;
endSamplePre = min(endSample, totalPre);
startSampleMain = max(startSample, totalPre + 1);
endSampleMain = min(endSample, totalMain + totalPre);

if ~onlyPulseInfoFlag
    % First, reproduce in each channel and silence the rest
    [x_pre, startPulseIndPre, endPulseIndPre] = successiveChannelSinusoids_inds( amplitude, phase, frequency, SampleRate, soundSamples, silenceSamples, startSamplePre, endSamplePre );
    
    % Then, reproduce everything at the same time
    t = (startSampleMain:endSampleMain)'/SampleRate;
    numSamplesMain = numel(t);
    x_main = zeros(numSamplesMain, numChann, numFreq);
    for c = 1:numChann
        for f = 1:numFreq
            x_main(:, c, f) = amplitude(c, f) * cos(2*pi*frequency(f)*t + phase(c, f));
        end
    end
    x_main = sum(x_main, 3);
    
    % Put both signals together
    x = [x_pre; x_main];
else
    x = [];
end

if numSamplesMain > 0
    startPulseIndMain = startSampleMain;
    endPulseIndMain = endSampleMain;
else
    startPulseIndMain = [];
    endPulseIndMain = [];
end

startInd = [startPulseIndPre; startPulseIndMain];
endInd = [endPulseIndPre; endPulseIndMain];

end