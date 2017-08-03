function [ pulseCoefMat, pulseLimits ] = coefficients2pulseSignalParameters( coefficients, frequency, SampleRate, option, varargin )
% coefficients. (numChannels x numFreq) matrix.
% frequency. numFreq-element vector.

numChann = size(coefficients, 1);
numFreq = size(coefficients, 2);

switch option
    case 'preludeAndMain'
       
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
        
        % Then, reproduce everything at the same time
        pulseCoefMat_main = permute(coefficients, [3, 1, 2]);
        pulseLimits = [pulseLimitsPre; [pulseLimitsPre(end) + silenceSamples, pulseLimitsPre(end) + silenceSamples + mainSamples]; ...
            [pulseLimitsPre(end) + silenceSamples + mainSamples, pulseLimitsPre(end) + silenceSamples + mainSamples + silenceSamples]];
        pulseCoefMat = [pulseCoefMat_pre; pulseCoefMat_main; zeros(1, numChann, numFreq)];
    case 'prelude'
        p = inputParser();
        
        addOptional(p, 'soundTime', 1);
        addOptional(p, 'silenceTime', 1);
        addOptional(p, 'numRepetitions', 1);
        
        parse(p, varargin{:});
        
        soundTime = p.Results.soundTime;
        silenceTime = p.Results.silenceTime;
        numRep = p.Results.numRepetitions;
        
        % Process       
        soundSamples = ceil(soundTime * SampleRate);
        silenceSamples = ceil(silenceTime * SampleRate);
        
        % Reproduce in each channel and silence the rest
        [pulseCoefMat, pulseLimits] = successiveChannelSinusoids( coefficients, frequency, soundSamples, silenceSamples, numRep);
    case 'main'
        % Prelude and main signal parameters
        p = inputParser();
        
        addOptional(p, 'soundTime', 1);
        addOptional(p, 'silenceTime', 1);
        addOptional(p, 'numRepetitions', 1);
        
        parse(p, varargin{:});
        
        soundTime = p.Results.soundTime;
        silenceTime = p.Results.silenceTime;
        numRep = p.Results.numRepetitions;  
        
        soundSamples = ceil(soundTime * SampleRate);
        silenceSamples = ceil(silenceTime * SampleRate);
        
        % Then, reproduce everything at the same time
        pulseCoefMat = repmat([permute(coefficients, [3, 1, 2]); zeros(1, numChann, numFreq)], numRep, 1);
        pulseLimits = repmat([0 soundSamples; soundSamples+1 soundSamples + silenceSamples], numRep + 1);
            
end




end