function H = OFDMchannelEstimation( transSignal, recSignal, coefMat, pulseLimits, freqs, sampleRate, varargin)

p = inputParser;
addParameter(p, 'margin', 0.5) % in seconds
parse(p, varargin{:})

margin = p.Results.margin;

% Check that each frequency is only transmitted by a channel simultaneously
aux = sum(coefMat ~= 0, 2);
assert(all(aux(:) <= 1), 'OFDMchannelEstimation:coefMatNotValid', 'The transmitted signal is not apt for OFDM channel estimation')

if ischar(transSignal)
    filename = transSignal;
    fromFile = true;
else
    fromFile = false;
end

[numPulses, numChannels, numFreqs] = size(coefMat);
numMicros = size(recSignal, 2);

freqsDif = diff(freqs);
if all(freqsDif == freqsDif(1))
    % The frequency steps are all equal
    freqStep = freqsDif(1);
else
    relFreqsDif = freqsDif/min(freqsDif);
    if all(relFreqsDif - round(relFreqsDif) == 0)
        % The frequency steps are all multiples of the minimum one
        freqStep = min(freqsDif);
    else
        error('OFDMchannelEstimation:wrongFreq', 'The vector of frequencies should have a constant frequency step')
    end
end
period = 1/freqStep;

transSpec = zeros(numChannels, numFreqs);
recSpec = zeros(numMicros, numChannels, numFreqs);
coefMatFlag = permute(coefMat ~= 0, [3 2 1]);
for p = 1:numPulses     
    startPulse = pulseLimits(p, 1);
    endPulse = pulseLimits(p, 2);
    
    durPulse = endPulse - startPulse;
    durFragment = period * floor((durPulse - margin*2)/period);
    durFragmentSamp = durFragment * sampleRate;
    
    startIndex = floor((startPulse + margin)*sampleRate) + 1;
    endIndex = startIndex + durFragmentSamp;
    
    % Calculate the transmitted signals spectra
    activeChannels = find(any(coefMat(p, :, :) ~= 0, 3)); % Active channels in this pulse
    if fromFile
        fragmentAllChann = audioread(filename, [startIndex, endIndex - 1], 'native');
    else
        fragmentAllChann = transSignal(startIndex:endIndex - 1, :);
    end
    for c = 1:length(activeChannels)
        chan = activeChannels(c);
        fragment = fragmentAllChann(:, chan);
        transSpec(chan, :) = freqz(fragment, 1, freqs, sampleRate);
    end
    
    % Calculate the received signals spectra    
    for m = 1:numMicros
        fragment = recSignal(startIndex:endIndex-1, m);
        X = freqz(fragment, 1, freqs, sampleRate); % (1 x numFreqs)
        
        for c = 1:length(activeChannels)
            chan = activeChannels(c);
            activeFreq = coefMatFlag(:, chan, p);
            recSpec(m, chan, activeFreq) = X(activeFreq);
        end
    end
    
    fprintf('%d/%d\n', p, numPulses);
end

flag = 0 ~= sum(coefMatFlag, 3)'; % This wouldn't be necessary for a full coefMat, but I expect the user to be a little retarded.
flagExt = repmat(permute(flag, [3 1 2]), [numMicros, 1, 1]);
transSpecExt = repmat(permute(transSpec, [3 1 2]), [numMicros, 1, 1]);

H = zeros(size(recSpec));
H(flagExt) = recSpec(flagExt)./transSpecExt(flagExt); % Channel frequency response

end