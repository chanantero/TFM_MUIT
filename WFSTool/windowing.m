function mask = windowing(type, varargin)

switch type
    case 'HanningModConstantRatio'
        p = inputParser;
        
        addOptional(p, 'numSamples', 1024)
        addOptional(p, 'constantRatio', 1)
        
        parse(p, varargin{:})
        
        numSamples = p.Results.numSamples;
        constantRatio = p.Results.constantRatio;
        
        mask = windowing1(numSamples, constantRatio);
        
    case 'HanningModRisingDuration'
        p = inputParser;
        
        addOptional(p, 'numSamples', 44100)
        addOptional(p, 'sampleRate', 44100)
        addOptional(p, 'risingDuration', 0.1)
        
        parse(p, varargin{:})
        
        numSamples = p.Results.numSamples;
        sampleRate = p.Results.sampleRate;
        risingDuration = p.Results.risingDuration;
        
        mask = windowing2(numSamples, sampleRate, risingDuration);
        
end


end

function mask = windowing1(numSamples, constantRatio)
hanningRatio = 1 - constantRatio;
numHann = floor(hanningRatio*numSamples);
numPre = ceil(numHann/2);
numPost = numHann - numPre;
hannWind = hann(numHann);

preInd = 1:numPre;
postInd = (numSamples-numPost+1):numSamples;
constantInd = numPre+1:numSamples-numPost;

mask = zeros(numSamples, 1);
mask(preInd) = hannWind(1:numPre);
mask(postInd) = hannWind(numPre+1:end);
mask(constantInd) = 1;

end

function mask = windowing2(numSamples, sampleRate, risingDuration)
risingSamples = floor(sampleRate*risingDuration);

numPre = min(risingSamples, floor(numSamples/2));
numPost = numPre;
hannWind = hann(numPre*2);

preInd = 1:numPre;
postInd = (numSamples-numPost+1):numSamples;
constantInd = numPre+1:numSamples-numPost;

mask = zeros(numSamples, 1);
mask(preInd) = hannWind(1:numPre);
mask(postInd) = hannWind(numPre+1:end);
mask(constantInd) = 1;

end