function [ pulseCoefMat, pulseLimits ] = coefficients2pulseSignalParameters( coefficients, frequency, SampleRate, option, varargin )
% coefficients. (numChannels x numFreq) matrix.
% frequency. numFreq-element vector.

numChann = size(coefficients, 1);
numFreq = size(coefficients, 2);

switch option
    case 'prelude'
       
        % Prelude and main signal parameters
        p = inputParser();
        
        addOptional(p, 'soundTime', 1);
        addOptional(p, 'silenceTime', 1);
        addOptional(p, 'numRepetitions', 5);
        addOptional(p, 'mainTime', 5);
        
        parse(p, varargin{:});
        
        soundTime = p.Results.soundTime;
        silenceTime = p.Results.silenceTime;
        numRep = p.Results.numRepetitions;
        mainTime = p.Results.mainTime;
        
        % Process       
        soundSamples = ceil(soundTime * SampleRate);
        silenceSamples = ceil(silenceTime * SampleRate);
        mainSamples = ceil(mainTime*SampleRate);
        
        % First, reproduce in each channel and silence the rest
        [pulseCoefMat_pre, pulseLimitsPre] = successiveChannelSinusoids( coefficients, frequency, soundSamples, silenceSamples, numRep);
        pulseLimitsPre = pulseLimitsPre - 1;
        
        % Then, reproduce everything at the same time
        pulseCoefMat_main = permute(coefficients, [3, 1, 2]);
        pulseLimits = [pulseLimitsPre; [pulseLimitsPre(end) + silenceSamples, pulseLimitsPre(end) + silenceSamples + mainSamples]; ...
            [pulseLimitsPre(end) + silenceSamples + mainSamples, pulseLimitsPre(end) + silenceSamples + mainSamples + silenceSamples]];
        pulseCoefMat = [pulseCoefMat_pre; pulseCoefMat_main; zeros(1, numChann, numFreq)];
        
end




end