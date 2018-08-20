%% Experiment 16.

% Created on 10/08/2018.

% Measures on the GTAC listening room. Reproduction and recording.

%% Preamble
pathSetUp;

imagesPath = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\TFM\Img\';

dataPathName = [globalPath, 'Data\'];
ID = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

%% System set up.
obj = SimulationController;
obj.ax.Parent.HandleVisibility = 'off';
obj.ax.HandleVisibility = 'off';

% Estimate visually the noise source position.
% Maket it reproduce with amplitude 1 and phase 0 in order to simplify
% things.
obj.NSposition = [3.35 -0.2 0]; % Assumed real position
obj.amplitude = 1;
obj.phase = 0;
obj.frequency = 440;

obj.WFSToolObj.setNumReceivers(1);

% As the audio driver only allows 96 channels and we need the reproduction
% of the noise source signal to be synchronized, it is necessary to
% dedicate one of the channels to the noise source. Of course, the corresponding
% array loudspeaker will be in silence.
% Set noise source channel through the GUI or use next line
NSchan = 1;
obj.WFSToolObj.noiseSourceChannelMapping(1) = NSchan;

%% Test that loudspeakers reproduce in an adequate way

% User change this
pulseDur = 1.5; % seconds
pulseSolap = 0.5; % seconds
testLaptop = true; % testLaptop should be false when reproducing in the listening room

% % % Don't touch the rest of the section % % %

% In order to test all loudspeakers work, we play the major musical scale
% back and forth (Do, Re, Mi, Fa, Sol, La, Si, Do', Si, La, Sol, Fa, Mi,
% Re), one note per loudspeaker, in an ordered way.
numWFS = 96;
if testLaptop
    coefMat = zeros(numWFS, 2, 8); 
    semiToneShift = [-9 -7 -5 -4 -2 0 2 3];
    freqLA = 440;
    freqs = freqLA .* 2.^(semiToneShift/12);
    prog = repmat([(1:7), (8:-1:2)], 1, ceil(numWFS/14));
    indFreq = prog(1:numWFS);
    indChan = mod(1:numWFS, 2)+1; 
    indPulse = 1:numWFS;
    inds = sub2ind([numWFS, 2, 8], indPulse, indChan, indFreq);
    coefMat(inds) = 1;
else % Test in GTAC listening room
    coefMat = zeros(numWFS, numWFS, 8); % zeros(numWFS, 2, 8);
    semiToneShift = [-9 -7 -5 -4 -2 0 2 3];
    freqLA = 440;
    freqs = freqLA .* 2.^(semiToneShift/12);
    prog = repmat([(1:7), (8:-1:2)], 1, ceil(numWFS/14));
    indFreq = prog(1:numWFS);
    indChan = 1:numWFS; % mod(1:numWFS, 2)+1;
    indPulse = 1:numWFS;
    inds = sub2ind([numWFS, numWFS, 8], indPulse, indChan, indFreq); % sub2ind([numWFS, 2, 8], indPulse, indChan, indFreq);
    coefMat(inds) = 1;
end

pulseStart = (0:numWFS-1)'*(pulseDur - pulseSolap);
pulseEnd = pulseStart + pulseDur;
pulseLimits = [pulseStart, pulseEnd];

sampleRate = 44100;
signalFunction = @(startSample, endSample) pulseCoefMat2signal(coefMat, pulseLimits,...
    freqs, sampleRate, startSample, endSample, 'type_pulseLimits', 'time');

% Debug
% ax = axes(figure);
% plot(ax, signalFunction(0, 44100*2))

obj.WFSToolObj.reproduceSignalFunction(signalFunction, sampleRate);

%% Generate frequency response signal

% In order to measure the acoustic paths, we are going to use the method of
% orthogonal frequecies.
% We are going to reproduce evenly spaced tones from 20Hz to 950Hz.
% However, every loudspeaker is going to reproduce different frequencies,
% so, in the reception of signals, we can, by spectral analysis,
% differenciate the signal from the different loudspeakers.

% Generate signals
numTotalChan = 96;
numSimultChan = 24; % Number of simultaneous channels reproducing
freqRange = [20, 950]; % Hz. It's a user preference, but it's not going to be exactly respectd
freqStep = 1; % Hz
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

% The time of the adquisition of a tone affect the signal to noise ratio?
% The longer the time, the higher the ratio. So, we can scale the tone
% coefficients in order to avoid saturation, and then compensate it by
% measuring during more time
maxPossible = repmat(sum(abs(coefMat), 3), [1, 1, numFreqs]);
coefMat(maxPossible ~= 0) = coefMat(maxPossible ~= 0)./maxPossible(maxPossible ~= 0);

pulseDur = (1/freqStep)*10 + 2;
silenceBetweenPulses = 3;
pulseStart = silenceBetweenPulses + (pulseDur + silenceBetweenPulses)*(0:numChanBlocks-1)';
pulseEnd = pulseStart + pulseDur;
pulseLimits = [pulseStart, pulseEnd];

signalStruct = struct('coefMat', coefMat, 'freqs', freqs, 'pulseLimits', pulseLimits, 'numSimultChan', numSimultChan, 'freqInd', freqInd);

sampleRate = 44100/4;
signalFunction = @(startSample, endSample) pulseCoefMat2signal(coefMat, pulseLimits,...
    freqs, sampleRate, startSample, endSample, 'type_pulseLimits', 'time');

% Write signal into file
audWriterObj = dsp.AudioFileWriter;
audWriterObj.Filename = [dataPathName, 'acousticPathReprodSignal_', ID, '.wav'];
audWriterObj.SampleRate = sampleRate;
audWriterObj.DataType = 'single';

% % Debug
% audReaderObj = dsp.AudioFileReader;
% filename = 'C:\Users\Rubén\Music\Salsa\Flor Pálida - Marc Anthony.mp3';
% audReaderObj.Filename = filename;
% audReaderObj.SamplesPerFrame = audReaderObj.SampleRate;
% 
% audWriterObj.SampleRate = audReaderObj.SampleRate;
% 
% count = 0;
% while ~isDone(audReaderObj)
%     signal = step(audReaderObj);
%     step(audWriterObj, single(signal));
%     count = count + 1;
% end
% 
% release(audReaderObj)
% release(audWriterObj);

sampIni = 1;
sampPerFrame = sampleRate;
sampEnd = sampIni + sampPerFrame - 1;
outOfRange = false;
count = 0; maxCount = ceil(pulseLimits(end)*sampleRate/sampPerFrame);
while outOfRange ~= 1
    [signal, outOfRange] = signalFunction(sampIni, sampEnd);
    sampIni = sampIni + sampPerFrame;
    sampEnd = sampEnd + sampPerFrame;
    signal(abs(signal) > 1) = sign(signal(abs(signal) > 1));
    step(audWriterObj, single(signal))
    count = count + 1;
    fprintf('count = %d/%d\n', count, maxCount);
end
release(audWriterObj);
save([dataPathName, 'signalStruct_', ID, '.mat'], '-struct', 'signalStruct');

% % Debug
ax = axes(figure);
signalEx = signalFunction(sampleRate*3, sampleRate*10);
plot(ax, signalEx(:, 1))
% [ax.Children(:).Visible] = deal('off');
% ax.Children(1).Visible = 'on';

%% Reproduce and record frequency response test signal
tic;
% Load frequency response test signal parameters
IDload = '2018-08-12_12-56-46';
signalStructLoaded = load([dataPathName, 'signalStruct_', IDload, '.mat']);
freqStep = signalStructLoaded.freqs(2) - signalStructLoaded.freqs(1);
T = 1/freqStep;
numChannels = 96;
numSimultChan = signalStructLoaded.numSimultChan;
pulseLimits = signalStructLoaded.pulseLimits;
freqs = signalStructLoaded.freqs;
numFreqs = length(freqs);
coefMat = signalStructLoaded.coefMat;
info = audioinfo([dataPathName, 'acousticPathReprodSignal_', IDload, '.wav']);
sampleRate = info.SampleRate;
filename = [dataPathName, 'acousticPathReprodSignal_', IDload, '.wav'];

simulatedReproduction = true;
if ~simulatedReproduction
    % Reproduce
    repRecObj = reproductorRecorder;
    repRecObj.setProps('mode', originType('file'), 1);
    repRecObj.setProps('enableProc', false);
    % filename = 'C:\Users\Rubén\Music\Salsa\Flor Pálida - Marc Anthony.mp3';
    repRecObj.setProps('audioFileName', filename, 1);
    repRecObj.executeOrder('play');
    recSignal = repRecObj.recorded{1}.';
else
    % Reproduce in simulation
    [signal, Fs] = audioread(filename, 'native');
    
    zPos = 1.65;
    WFSarrayOffset = [0.46 2.21 zPos]; % [x, y, z] coordinates. Useful for generating acoustic path IR.
    roomDim = [4.48, 9.13, 2.64];
    fs = sampleRate;
    c = 340;
    
    extRectXmin = min(obj.WFSposition(:, 1));
    extRectXmax = max(obj.WFSposition(:, 1));
    extRectYmin = min(obj.WFSposition(:, 2));
    extRectYmax = max(obj.WFSposition(:, 2));
    centerX = (extRectXmax + extRectXmin)/2;
    centerY = (extRectYmax + extRectYmin)/2;
    recPositions = [centerX, centerY, 0; centerX, centerY + 1, 0];
    
    % Room characteristics and impulse response of chamber
    beta = 0; % Average reflection coefficient of the walls of the chamber
    WFS_AcPath_previously_calculated = true;
    NS_AcPath_previously_calculated = true;
    appendFreeSpaceAcPaths = false;
    
    % Irrelevant variables that are necessary in order to SetupParameterScript
    % to work
    amplitude = 1;
    phase = 0;
    freqFilters = {};
    NSpositions = obj.NSRposition;
    
    SetupParametersScript
    AcousticPathCalculationScript
    
    obj.domain = 'time';
    
    WFSacPathIR = permute(WFS_IR(:, :, :, 1), [1, 3, 2]);
    NSacPathIR = repmat(permute(NS_IR, [1, 3, 2]), [1 2 1]);
    obj.setAcousticPaths('WFS', WFSacPathIR);
    obj.setAcousticPaths('NS', NSacPathIR);
    
    wfsChan = obj.WFSToolObj.loudspeakerMapping(1).destinationInd; % Not pretty sure
    nsChan = obj.WFSToolObj.loudspeakerMapping(2).destinationInd;
    fragLength = 2^12; % fragment length in samples
    numSamp = size(signal, 1);
    numFrag = ceil(numSamp/fragLength); % number of fragments
    previousSufix = zeros(numMicro, numSampIR - 1);
    recSignal = zeros(numMicro, numSamp); % Recorded signal. (numMicros x numSamples)
    for f = 1:numFrag
        fprintf('%d/%d\n', f, numFrag);
        selSamp = 1 + fragLength*(f - 1):min(fragLength*f, numSamp); % The min() is for the last iteration
        frag = [signal(selSamp, :)', zeros(numChannels, numSampIR-1)];
        obj.WFSToolObj.WFSarrayCoefficient = zeros(obj.numWFS, size(frag, 2));
        obj.WFSToolObj.WFSarrayCoefficient(obj.WFSToolObj.loudspeakerMapping(1).originInd, :) = frag(wfsChan, :);
        obj.WFSToolObj.noiseSourceCoefficient_complete = zeros(2, size(frag, 2));
        obj.WFSToolObj.noiseSourceCoefficient_complete(obj.WFSToolObj.loudspeakerMapping(2).originInd, :) = frag(nsChan, :);
        obj.WFSToolObj.simulTheo.updateField();
        recFrag = obj.WFSToolObj.simulField;
        recFrag(:, 1:numSampIR-1) = recFrag(:, 1:numSampIR-1) + previousSufix; % Add previous sufix
        previousSufix = recFrag(:, end - (numSampIR - 1) + 1:end); % Set sufix for next loop iteration
        recSignal(:, selSamp) = recFrag(:, 1:length(selSamp)); % Save the signal without the sufix. We use length(selSamp) and not fragLength for the last iteration
    end
end

% Crop signal with an adequate duration
numMicros = size(recSignal, 1);
margin = 1;

numPulses = size(pulseLimits, 1);
transSpec = zeros(numChannels, numFreqs);
recSpec = zeros(numMicros, numChannels, numFreqs);
coefMatFlag = permute(coefMat ~= 0, [3 2 1]);
for p = 1:numPulses
    start = pulseLimits(p, 1) + margin;
    ending = start + floor((pulseLimits(p, 2) - margin - start)/T);
    
    startIndex = floor(start*sampleRate) + 1;
    endIndex = floor(ending*sampleRate) + 1;
    
    % Calculate the transmitted signals spectra
    activeChannels = find(any(coefMat(p, :, :) ~= 0, 3)); % Active channels in this pulse
    for c = 1:length(activeChannels)
        fragment = signal(startIndex:endIndex, activeChannels(c));
        transSpec(activeChannels(c), :) = freqz(fragment, 1, freqs, sampleRate);
    end
    
    % Calculate the received signals spectra    
    for m = 1:numMicros
        fragment = recSignal(m, startIndex:endIndex);
        X = freqz(fragment, 1, freqs, sampleRate); % (1 x numFreqs)
        
        for c = 1:numChannels
            activeFreq = coefMatFlag(:, c, p);
            recSpec(m, c, activeFreq) = X(activeFreq);
        end       
    end
end
FR = recSpec./repmat(permute(transSpec, [3 1 2]), [numMicros, 1, 1]);

ax = axes(figure);
plot(ax, freqs, squeeze(abs(FR(1, 1:2, :))))
plot(ax, freqs, squeeze(abs(recSpec(1, 1:2, :))))
plot(ax, freqs, squeeze(abs(transSpec(1:2, :))))

% Check the result is equal to the old form of doing it. All must be
% simulated (by now) of course.
if simulatedReproduction
    wfsInd = obj.WFSToolObj.loudspeakerMapping(1).originInd;
    wfsToLoudsInd = obj.WFSToolObj.loudspeakerMapping(1).destinationInd;
    nsInd = obj.WFSToolObj.loudspeakerMapping(2).originInd;
    nsToLoudsInd = obj.WFSToolObj.loudspeakerMapping(2).destinationInd;
    FRtheo = zeros(numMicros, numChannels, numFreqs);
    FRtheo(:, wfsToLoudsInd, :) = permute(WFS_FR(:, :, wfsInd), [1, 3, 2]);
    FRtheo(:, nsToLoudsInd, :) = permute(NS_FR(:, :, nsInd), [1, 3, 2]);
    recFlag = repmat(any(coefMat ~= 0, 1), [numMicros, 1, 1]);
    rel = ones(size(FR));
    rel(recFlag) = FR(recFlag)./FRtheo(recFlag);
    fig = figure;
    ax = axes(fig);
    histogram(ax, abs(rel(recFlag)))
    histogram(angle(rel(recFlag)))
    max(abs(rel(recFlag)))
    min(abs(rel(recFlag)))
    max(angle(rel(recFlag)))
    min(angle(rel(recFlag)))
    % And they are!!! :D
    WFS_FR_a = WFS_FR;
    NS_FR_a = NS_FR;
end

% Each channel has associated a vector of frequencies and a vector of
% frequency responses
% We can complete missing frequencies by interpolating
FRint = zeros(size(recSpec));
for m = 1:numMicros
    for c = 1:numChannels
        nonZero = recSpec(m, c, :) ~= 0;
        magn = interp1(freqs(nonZero), abs(squeeze(FR(m, c, nonZero))), freqs, 'linear');
        phase = interp1(freqs(nonZero), unwrap(angle(squeeze(FR(m, c, nonZero)))), freqs, 'linear');
        FRint(m, c, :) = magn.*exp(1i*phase);
    end
end
toc1 = toc;
%% Estimation of results based on measured acoustic path responses

% ---- No optimization ----
% Set parameters
NSpositions = [3.35 -0.2 0]; % Assumed real position
amplitude = 1;
phase = 0;
freqFilters = {};
freqFiltDelays = [];

simulScriptFlag = false;
if simulScriptFlag
    fs = sampleRate;
    c = 340;
    
    extRectXmin = min(obj.WFSposition(:, 1));
    extRectXmax = max(obj.WFSposition(:, 1));
    extRectYmin = min(obj.WFSposition(:, 2));
    extRectYmax = max(obj.WFSposition(:, 2));
    centerX = (extRectXmax + extRectXmin)/2;
    centerY = (extRectYmax + extRectYmin)/2;
    recPositions = [centerX, centerY, 0; centerX, centerY + 1, 0]; % The actual positions don't really matter, but the number of receiver positions must be the same as microphones
    
    % Room characteristics and impulse response of chamber
    beta = 0;
    WFS_AcPath_previously_calculated = true;
    NS_AcPath_previously_calculated = true;
    wfsChan = obj.WFSToolObj.loudspeakerMapping(1).destinationInd; % Not pretty sure
    nsChan = obj.WFSToolObj.loudspeakerMapping(2).destinationInd;
    WFS_FR = permute(FRint, [1, 3, 2]);
    WFS_FR(:, :, nsChan) = 0;
    NS_FR = permute(FRint(:, nsChan, :), [1, 3, 2]);

    % WFS options
    frequencyCorrection = true;
    attenuationType = 'Ruben';
    
    % Simulation options
    timeDomainActive = false;
    fakeTimeProcessing = false;
    frequencyDomainActive = true;
    automaticLengthModification = false;
    
    SetupParametersScript
    simulationScript
    
    % recNScoef_freq
    % recCoef_freq
    % recWFScoef_freq
    % WFScoef_freq
        
    sSimul = s';

    [sExtTheo, corrFactInd, corrFactGlob, gainInd, gainGlob, corrFactAver, gainAver] =...
        SimulationController.addCancellationParametersToStructure(sSimul);
    
    ax = axes(figure);
    plot(ax, freqs, 10*log10(gainAver));
    plot(ax, freqs, abs(corrFactAver))
    
else
    % No completado
    wfsChan = obj.WFSToolObj.loudspeakerMapping(1).destinationInd; % Not pretty sure
    nsChan = obj.WFSToolObj.loudspeakerMapping(2).destinationInd;
    WFS_FR = FRint;
    WFS_FR(:, nsChan, :) = 0;
    NS_FR = FRint(:, nsChan, :);
    
    % WFS calculation of WFS coefficients and simulation
    obj.NSposition = NSpositions;
    obj.domain = 'frequency';
    fieldNS = zeros(numMicros, numFreqs);
    fieldWFS = zeros(numMicros, numFreqs);
    WFScoefs = zeros(obj.numWFS, numFreqs);
    
    for f = 1:numFreqs
        obj.frequency = freqs(f);
        obj.WFSToolObj.WFScalculation();
        WFScoef = obj.WFScoef;
        WFScoefs(:, f) = WFScoef;
        fieldNS(:, f) = NS_FR(:, :, f)*obj.NSRcoef;
        fieldWFS(:, f) = WFS_FR(:,:,f)*WFScoef;
    end
    field = fieldNS + fieldWFS;
        
    % Calculation of gain
    % Make it automatic by generating an s structure and using
    s = repmat(obj.generateBasicExportStructure, numFreqs, 1);
    for f = 1:numFreqs
        s(f).recCoef = field(:, f);
        s(f).recNScoef = fieldNS(:, f);
        s(f).recWFScoef = fieldWFS(:, f);
        s(f).WFScoef = WFScoefs(:, f);
    end
    [sExt, corrFactInd, corrFactGlob, gainInd, gainGlob, corrFactAver, gainAver] =...
        SimulationController.addCancellationParametersToStructure(s);
    
    sSimul = s;
    
end

% ---- Volume optimization ----
% Choose frequencies high enough so we are not in the low-frequency zone
minFreq = 500;
maxFreq = 850;

recNScoef = [sSimul.recNScoef];
recWFScoef = [sSimul.recWFScoef];

sel = freqs >= minFreq & freqs <= maxFreq;
fieldWFSaux = recWFScoef(:, sel);
fieldNSaux = recNScoef(:, sel);
corrFact = -fieldWFSaux(:)\fieldNSaux(:);
corrVol = real(corrFact);

sSimulCorrVol = sSimul;
for f = 1:numFreqs
    sSimulCorrVol(f).WFScoef = sSimulCorrVol(f).WFScoef*corrVol;
    sSimulCorrVol(f).recWFScoef = sSimulCorrVol(f).recWFScoef*corrVol;
    sSimulCorrVol(f).recCoef = sSimulCorrVol(f).recWFScoef + sSimulCorrVol(f).recNScoef;
end

% ---- Volume and noise source location optimization ----
% Completar
% selFreqs = find(sel);
% for f = 1:length(selFreqs)
%     obj.frequency = selFreqs(f);
%     
%     obj.findBestVirtualSourceParameters2;
% end

%% Reproduction and recording of different cases

% Define chirp signal
durSign = 4; % Duration of tone for time processing
t = (0:ceil(durSign*sampleRate)-1)/sampleRate;
NSsignal = chirp(t, min(freqs), durSign, max(freqs) + 50);
prefixDuration = 2;
prefixNumSamples = floor(prefixDuration*sampleRate) + 1;
NSsignal = [zeros(1, prefixNumSamples), NSsignal, zeros(1, prefixNumSamples )];
t = (0:length(NSsignal)-1)/sampleRate;

% Frequency filters
magnFiltOrder = 2^11;
hilbertFiltOrder = 2^13;
[freqFilter, freqFiltDelay] = getFrequencyFilter(magnFiltOrder, hilbertFiltOrder, sampleRate);

% Calculate WFS signals
obj.domain = 'time';
obj.NScoef = NSsignal;
obj.NSVcoef = -NSsignal;

obj.WFSToolObj.freqFilter = freqFilter;

freqFilterResp = freqz(freqFilter, 1, freqs, sampleRate);
ax = axes(figure);
plot(ax, freqs, abs(freqFilterResp), freqs, sqrt(freqs/340))
plot(ax, freqs, rad2deg(wrapToPi(angle(freqFilterResp) + freqFiltDelay/sampleRate*2*pi*freqs)))

obj.WFSToolObj.frequencyCorrection = true;
obj.WFSToolObj.attenuationType = 'Ruben'; attenuationType = 'Ruben';
obj.WFSToolObj.updateFiltersWFS();
obj.WFSToolObj.automaticLengthModification = false;
obj.WFSToolObj.WFScalculation();

numSamp = length(NSsignal);
customSignal = zeros(numSamp, numChannels);
wfsIndOrig = obj.WFSToolObj.loudspeakerMapping(1).originInd;
wfsIndDest = obj.WFSToolObj.loudspeakerMapping(1).destinationInd;
nsIndDest = obj.WFSToolObj.loudspeakerMapping(2).destinationInd;
customSignal(:, nsIndDest) = NSsignal(:);
customSignal(:, wfsIndDest) = obj.WFScoef(wfsIndOrig, :).';

% Reproduce and record:
if ~simulatedReproduction
    repRecObj = reproductorRecorder;
    repRecObj.setProps('mode', originType('custom'), 1);
    repRecObj.setProps('enableProc', false);
    repRecObj.setProps('Fs_player', sampleRate, 1);
    repRecObj.setProps('Fs_recorder', sampleRate, 1);
    
% ---- No optimization ----
    
    % - Only noise source
    aux = customSignal;
    aux(:, wfsIndDest) = 0;
    repRecObj.setProps('customSignal', aux, 1);
    repRecObj.executeOrder('play');
    recSignalNS = repRecObj.recorded{1}.';
    
    % - Only WFS
    aux = customSignal;
    aux(:, nsIndDest) = 0;
    repRecObj.setProps('customSignal', aux);
    repRecObj.executeOrder('play');
    recSignalWFS = repRecObj.recorded{1}.';
    
    % - All
    repRecObj.setProps('customSignal', customSignal);
    repRecObj.executeOrder('play');
    recSignal = repRecObj.recorded{1}.';
    
% ---- Optimized volume ----
    customSignal(:, wfsIndDest) = corrVol*obj.WFScoef(wfsIndOrig, :).';
    
    % - Only noise source
    aux = customSignal;
    aux(:, wfsIndDest) = 0;
    repRecObj.setProps('customSignal', aux);
    repRecObj.executeOrder('play');
    recSignalNScorrVol = repRecObj.recorded{1}.';
    
    % - Only WFS
    aux = customSignal;
    aux(:, nsIndDest) = 0;
    repRecObj.setProps('customSignal', aux);
    repRecObj.executeOrder('play');
    recSignalWFScorrVol = repRecObj.recorded{1}.';
    
    % - All
    repRecObj.setProps('customSignal', customSignal);
    repRecObj.executeOrder('play');
    recSignalCorrVol = repRecObj.recorded{1}.';
    
% ---- Optimized volume and noise source location ----
    % No completado
else
    
    % Room characteristics and impulse response of chamber
    beta = 0; % Average reflection coefficient of the walls of the chamber
    WFS_AcPath_previously_calculated = true;
    NS_AcPath_previously_calculated = true;
    appendFreeSpaceAcPaths = false;
    automaticLengthModification = false;
        
    SetupParametersScript
    AcousticPathCalculationScript
    
    obj.domain = 'time';
    
    WFSacPathIR = permute(WFS_IR(:, :, :, 1), [1, 3, 2]);
    NSacPathIR = repmat(permute(NS_IR, [1, 3, 2]), [1 2 1]);
    obj.setAcousticPaths('WFS', WFSacPathIR);
    obj.setAcousticPaths('NS', NSacPathIR);
    
% ---- No optimization ----
    
    % - Only noise source
    obj.WFSToolObj.real = [true; false];
    obj.WFSToolObj.virtual = [false; false];
    obj.WFSToolObj.WFScalculation();
    obj.WFSToolObj.simulate;
    recSignalNS = obj.microCoef;
    obj.WFSToolObj.virtual = [false; true];
    
    % - Only WFS
    obj.WFSToolObj.real = [false; false];
    obj.WFSToolObj.WFScalculation();
    obj.WFSToolObj.simulate;
    recSignalWFS = obj.microCoef;
    obj.WFSToolObj.real = [true; false];
    
    % Debug
    ax = axes(figure);
    plot(ax, t, NSsignal, t, obj.WFScoef(89, :))
    nsSignalFreq = freqz(NSsignal, 1, freqs, sampleRate);
    wfsSignalFreq = freqz(obj.WFScoef(89, :), 1, freqs, sampleRate);
    ax = axes(figure);
    dist = calcDistances(obj.NSRposition, obj.WFSposition(89,:));
    cosAlpha = dot((obj.WFSposition(89,:) - obj.NSRposition), [0 1 0]);
    plot(ax, freqs, abs(wfsSignalFreq./nsSignalFreq))
    plot(ax, freqs, rad2deg(wrapToPi(angle(wfsSignalFreq./nsSignalFreq) + dist/340*2*pi*freqs)))
    plot(ax, freqs, angle(wfsSignalFreq./nsSignalFreq) )
    
    % size(obj.filtersWFS_IR) (numWFS x numNS x numSamp)
    filt = squeeze(obj.WFSToolObj.filtersWFS_IR(89, 2, :));
    plot(ax, filt)
    filtSpec = freqz(filt, 1, freqs, sampleRate);
    plot(ax, freqs, angle(filtSpec))
    
    
    % - All
    obj.WFSToolObj.WFScalculation();
    obj.WFSToolObj.simulate;
    recSignal = obj.microCoef;

% ---- Optimized volume ----
    
    % - Only noise source
    obj.WFSToolObj.virtual = [false; false];
    obj.WFSToolObj.WFScalculation();
    obj.WFSToolObj.simulate;
    recSignalNScorrVol = obj.microCoef;
    obj.WFSToolObj.virtual = [false; true];
    
    % - Only WFS
    obj.WFSToolObj.real = [false; false];
    obj.WFSToolObj.WFScalculation();
    obj.WFScoef = obj.WFScoef*corrVol;
    obj.WFSToolObj.simulate;
    recSignalWFScorrVol = obj.microCoef;
    obj.WFSToolObj.real = [true; false];
    
    % - All
    obj.WFSToolObj.WFScalculation();
    obj.WFScoef = obj.WFScoef*corrVol;
    obj.WFSToolObj.simulate;
    recSignalCorrVol = obj.microCoef;
    
end



%% Comparison of measures and estimations. They should be similar.

% Get the received signal frequency spectrum
oper = @(x) freqz(x, 1, freqs, sampleRate);

recNS = oneDimOperOverMultiDimArray(oper, recSignalNS, 2);
recWFS = oneDimOperOverMultiDimArray(oper, recSignalWFS, 2);
rec = oneDimOperOverMultiDimArray(oper, recSignal, 2);
recNScorrVol = oneDimOperOverMultiDimArray(oper, recSignalNScorrVol , 2);
recWFScorrVol = oneDimOperOverMultiDimArray(oper, recSignalWFScorrVol , 2);
recCorrVol  = oneDimOperOverMultiDimArray(oper, recSignalCorrVol , 2);

obj.domain = 'frequency';
sExp = repmat(obj.generateBasicExportStructure, numFreqs, 1);
sExpCorrVol = repmat(obj.generateBasicExportStructure, numFreqs, 1);
for f = 1:numFreqs
    sExp(f).recCoef = rec(:, f);
    sExp(f).recNScoef = recNS(:, f);
    sExp(f).recWFScoef = recWFS(:, f);
    sExp(f).Frequency = freqs(f);
    
    sExpCorrVol(f).recCoef = recCorrVol(:, f);
    sExpCorrVol(f).recNScoef = recNScorrVol(:, f);
    sExpCorrVol(f).recWFScoef = recWFScorrVol(:, f);
    sExpCorrVol(f).Frequency = freqs(f);
end
[sExt, corrFactInd, corrFactGlob, gainInd, gainGlob, corrFactAver, gainAver] =...
        SimulationController.addCancellationParametersToStructure(sExp); 
[sExtCorrVol, corrFactIndCorrVol, corrFactGlobCorrVol, gainIndCorrVol, gainGlobCorrVol, corrFactAverCorrVol, gainAverCorrVol] =...
        SimulationController.addCancellationParametersToStructure(sExpCorrVol);
    
% Compare it with the simulated one according to the measured acoustic
% paths
[sExtSimul, corrFactIndSimul, corrFactGlobSimul, gainIndSimul, gainGlobSimul, corrFactAverSimul, gainAverSimul] =...
        SimulationController.addCancellationParametersToStructure(sSimul);
[sExtSimulCorrVol, corrFactIndSimulCorrVol, corrFactGlobSimulCorrVol, gainIndSimulCorrVol, gainGlobSimulCorrVol, corrFactAverSimulCorrVol, gainAverSimulCorrVol] =...
        SimulationController.addCancellationParametersToStructure(sSimulCorrVol);

ax = axes(figure, 'NextPlot', 'Add');
plot(ax, freqs, 10*log10(gainAver))
plot(ax, freqs, 10*log10(gainAverSimul))
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = 'Average gain (dB)';

ax = axes(figure, 'NextPlot', 'Add');
plot(ax, freqs, 10*log10(gainAverCorrVol))
plot(ax, freqs, 10*log10(gainAverSimulCorrVol))
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = 'Average gain (dB)';

% Debug
% The received signal spectrum from the noise source is the same
NSspec = freqz(NSsignal, 1, freqs, sampleRate);
resp = recNS./repmat(NSspec, 2, 1);
recNSsimul = [sSimul.recNScoef];
ax = axes(figure);
plot(ax, freqs, abs(resp), freqs, abs(recNSsimul))

% The received signal spectrum from the WFS array is the same? No, but it's
% similar. I think the differece is due to the fact that FRint is not
% exactly as the simulated acoustic paths, since it is an interpolation. Or
% maybe it is because the frequency response of the frequency filter is not
% ideal
resp = recWFS./repmat(NSspec, 2, 1);
recNSsimul = [sSimul.recWFScoef];
ax = axes(figure);
plot(ax, freqs, abs(resp), freqs, abs(recNSsimul)) % plot(ax, freqs, angle(resp), freqs, angle(recNSsimul))

ax = axes(figure);
plot(ax, t, recSignalNS(1,:), t, recSignalWFS(1,:))

ax

%% View of results

% The transmitted signal by the noise source is:
    % Time representation
    ax = axes(figure);
    plot(ax, t, NSsignal)
    ax.XLabel.String = 'Time (s)';
    ax.YLabel.String = '$\signal[ns][time]$'; ax.YLabel.Interpreter = 'latex';
    Plot2LaTeX(ax.Parent, [imagesPath, 'Experiment16_NSsignalTime'])

    % Frequency representation
    ax = axes(figure);
    plot(ax, freqs, abs(NSspec));
    ax.XLabel.String = 'Frequency (Hz)';
    ax.YLabel.String = '$\signal[ns][frequency]$'; ax.YLabel.Interpreter = 'latex';
    
% The received signals from the noise source in both microphones are:
    % Time representation
    % Frequency representation
    
% The received total signals are:
    % ---- No optimization ----
        % Time representation
        % Frequency representation
    
    % ---- Volume optimization ----    
        % Time representation
        % Frequency representation
        
        
        