function varargout = OFDMchannelEstimationSignal(numTotalChan, numSimultChan, freqRange, freqStep, pulseDurationInPeriods, silenceBetweenPulses, varargin)

p = inputParser;
addOptional(p, 'sampleRate', 44100)
addParameter(p, 'filename', [])
parse(p, varargin{:})

filename = p.Results.filename;
sampleRate = p.Results.sampleRate;

writeToFile = ~isempty(filename);

% Generate signals
freqStepChan = freqStep*numSimultChan; % Frequency step of each channel
numBlocks = ceil((freqRange(2) - freqRange(1))/freqStepChan);
freqs = (0:numBlocks*numSimultChan-1)*freqStep + freqRange(1);
freqInd = reshape(1:numBlocks*numSimultChan, [numSimultChan, numBlocks])';
numFreqs = numBlocks*numSimultChan;

numChanBlocks = ceil(numTotalChan/numSimultChan);
coefMat = zeros(numChanBlocks, numTotalChan, numFreqs);
for chanBlock = 1:numChanBlocks
    simultChanInd = (1:numSimultChan) + numSimultChan*(chanBlock - 1);
    simultChanInd(simultChanInd > numTotalChan) = [];
    subindSimultChan = repmat(simultChanInd, [numBlocks, 1]);

    inds = sub2ind([numChanBlocks, numTotalChan, numFreqs], chanBlock*ones(size(subindSimultChan)), subindSimultChan, freqInd(:, 1:length(simultChanInd)));
    coefMat(inds) = exp(1i*rand(size(inds))*2*pi);
end

maxPossible = repmat(sum(abs(coefMat), 3), [1, 1, numFreqs]);
coefMat(maxPossible ~= 0) = 0.5*coefMat(maxPossible ~= 0)./maxPossible(maxPossible ~= 0);

prefixSufixDuration = 1;
period = 1/freqStep;
pulseDur = period*pulseDurationInPeriods + prefixSufixDuration*2;
pulseStart = silenceBetweenPulses + (pulseDur + silenceBetweenPulses)*(0:numChanBlocks-1)';
pulseEnd = pulseStart + pulseDur;
pulseLimits = [pulseStart, pulseEnd];

signalStruct = struct('coefMat', coefMat, 'freqs', freqs, 'pulseLimits', pulseLimits, 'freqInd', freqInd);
if writeToFile
    save([filename, '.mat'], '-struct', 'signalStruct');

    % Write signal into file
    signalFunction = @(startSample, endSample) pulseCoefMat2signal(coefMat, pulseLimits,...
        freqs, sampleRate, startSample, endSample, 'type_pulseLimits', 'time');
    
    audWriterObj = dsp.AudioFileWriter;
    audWriterObj.Filename = [filename, '.wav'];
    audWriterObj.SampleRate = sampleRate;
    audWriterObj.DataType = 'single';
    
    sampIni = 1;
    sampPerFrame = sampleRate;
    sampEnd = sampIni + sampPerFrame - 1;
    outOfRange = false;
    % count = 0; maxCount = ceil(pulseLimits(end)*sampleRate/sampPerFrame);
    while outOfRange ~= 1
        [signal, outOfRange] = signalFunction(sampIni, sampEnd);
        sampIni = sampIni + sampPerFrame;
        sampEnd = sampEnd + sampPerFrame;
        signal(abs(signal) > 1) = sign(signal(abs(signal) > 1));
        step(audWriterObj, single(signal))
        %     count = count + 1;
        %     fprintf('count = %d/%d\n', count, maxCount);
    end
    release(audWriterObj);
    
    if nargout == 1
        varargout = {[]};
    else
        varargout = {};
    end
else
    signal = pulseCoefMat2signal(coefMat, pulseLimits, freqs, sampleRate, 0, pulseLimits(end, 2), 'type_pulseLimits', 'time', 'type_marker', 'time');
    varargout = {signal, signalStruct};
end

end

