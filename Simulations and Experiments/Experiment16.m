%% Experiment 16.

% Created on 10/08/2018.

% Measures on the GTAC listening room. Reproduction and recording.

%% Preamble
pathSetUp;

imagesPath = 'C:\Users\Rub�n\Google Drive\Telecomunicaci�n\M�ster 2� Curso 2015-2016\TFM MUIT\Documentos\Img\';

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
% filename = 'C:\Users\Rub�n\Music\Salsa\Flor P�lida - Marc Anthony.mp3';
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
    % filename = 'C:\Users\Rub�n\Music\Salsa\Flor P�lida - Marc Anthony.mp3';
    repRecObj.setProps('audioFileName', filename, 1);
    repRecObj.executeOrder('play');
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
    recFlag = any(coefMat ~= 0, 1);
    rel = ones(size(FR));
    rel(recFlag) = FR(recFlag)./FRtheo(recFlag);
    histogram(abs(rel(recFlag)))
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

%% Estimation of results based on measured acoustic path responses

% ---- No optimization ----
% Set parameters
NSpositions = [3.35 -0.2 0]; % Assumed real position
amplitude = 1;
phase = 0;
freqFilters = {};
freqFiltDelays = [];

simulScriptFlag = true;
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
    
    [sExt, corrFactInd, corrFactGlob, gainInd, gainGlob, corrFactAver, gainAver] =...
        SimulationController.addCancellationParametersToStructure(s);
    
%     ax = axes(figure);
%     plot(ax, freqs, 10*log10(gainAver));
%     plot(ax, freqs, abs(corrFactAver))
    
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
    field = zeros(1, numFreqs);
    fieldNS = zeros(1, numFreqs);
    fieldWFS = zeros(1, numFreqs);
    
    for f = 1:numFreqs
        obj.frequency = freqs(f);
        obj.WFSToolObj.WFScalculation();
        WFScoef = obj.WFScoef;
        fieldNS(f) = NS_FR(:, :, f)*obj.NSRcoef;
        fieldWFS(f) = WFS_FR(:,:,f)*WFScoef;
        field(f) = fieldNS(f) + fieldWFS(f);
    end
        
    % Calculation of gain
    % Make it automatic by generating an s structure and using
    s = repmat(obj.generateBasicExportStructure, numFreqs, 1);
    for f = 1:numFreqs
        s(f).recCoef = field(f);
        s(f).recNScoef = fieldNS(f);
        s(f).recWFScoef = fieldWFS(f);
    end
    [sExt, corrFactInd, corrFactGlob, gainInd, gainGlob, corrFactAver, gainAver] =...
        SimulationController.addCancellationParametersToStructure(s);   
end

% ---- Volume optimization ----
% Choose frequencies high enough so we are not in the low-frequency zone
minFreq = 500;
maxFreq = 850;

sel = freqs >= minFreq & freqs <= maxFreq;
fieldWFSaux = recWFScoef_freq(:, sel);
fieldNSaux = recNScoef_freq(:, sel);
corrFact = -fieldWFSaux(:)\fieldNSaux(:);
corrVol = real(corrFact);

% ---- Volume and noise source location optimization ----
% Completar
% selFreqs = find(sel);
% for f = 1:length(selFreqs)
%     obj.frequency = selFreqs(f);
%     
%     obj.findBestVirtualSourceParameters2;
% end

%% Reproduction and recording of different cases
% s should be of size (numMicros, numFreqs)

% ---- No optimization ----

% Define chirp signal
freqs = 0;
durSign = 4; % Duration of tone for time processing
t = (0:ceil(durSign*fs)-1)/fs;
NSsignal = chirp(t, min(freqs), durSign, max(freqs));
% Simulate only the noise source
obj.WFSToolObj.virtual = [false; false];
obj.WFSToolObj.freqFilter = 1;
obj.WFSToolObj.WFScalculation();
obj.WFSToolObj.simulate();
recNS_signal = obj.WFSToolObj.simulField;
recNScoef_time(:, f, ns, rt) = exp(1i*preDelayPhaseShift) * signal2pulseCoefficientMatrix([0 durSign], fcurr, 1, recNS_signal', fs).';

recNS_signals(:, :, f, ns, rt) = recNS_signal;
                    
obj.WFSToolObj.virtual = [false; true];
% Simulate all together
obj.WFSToolObj.freqFilter = freqFilters{filt};
obj.WFSToolObj.updateFiltersWFS();
obj.WFSToolObj.WFScalculation();
obj.WFSToolObj.simulate();
rec_signal = obj.WFSToolObj.simulField;

% Identify IQ component
recWFScoef_time(:, f, ns, rt, filt) = exp(1i*preDelayPhaseShift) * signal2pulseCoefficientMatrix([0 durSign], fcurr, 1, (rec_signal - recNS_signal)', fs);
recCoef_time(:, f, ns, rt, filt) = exp(1i*preDelayPhaseShift) * signal2pulseCoefficientMatrix([0 durSign], fcurr, 1, rec_signal', fs);
WFScoef_time(:, f, ns, rt, filt) = exp(1i*preDelayPhaseShift) * signal2pulseCoefficientMatrix([0 durSign], fcurr, 1, obj.WFSToolObj.WFSarrayCoefficient', fs);

rec_signals(:, :, f, ns, rt, filt) = rec_signal;
                        

% 3 pulses:
% - Only noise.
% - Only WFS.
% - Both.
coefMat = zeros(3, numChannels, numFreqs); % (3 x 96 x numFreqs)
WFScoefMat = [s.WFScoef];
NScoefMat = [s.NSVcoef];
nsIndOrig = obj.WFSToolObj.loudspeakerMapping(2).originInd;
wfsIndOrig = obj.WFSToolObj.loudspeakerMapping(1).originInd;
coefMat(1, nsChan, :) = permute(NScoefMat(nsIndOrig, :), [3 1 2]);
coefMat(2, wfsChan, :) = permute(WFScoefMat(wfsIndOrig, :), [3 1 2]);
coefMat(3, nsChan, :) = permute(NScoefMat(nsIndOrig, :), [3 1 2]);
coefMat(3, wfsChan, :) = permute(WFScoefMat(wfsIndOrig, :), [3 1 2]);

pulseLimits = [1 3; 4 6; 7 9];

sampleRate = 44100/4;
signalFunction = @(startSample, endSample) pulseCoefMat2signal(coefMat, pulseLimits,...
    freqs, sampleRate, startSample, endSample, 'type_pulseLimits', 'time');

obj.WFSToolObj.reproduceSignalFunction(signalFunction, sampleRate);
recField = obj.WFSToolObj.recorded;

% ---- Optimized volume ----
coefMat = zeros(3, numChannels, numFreqs); % (3 x 96 x numFreqs)
WFScoefMat = [s.WFScoef];
NScoefMat = [s.NSVcoef];
coefMat(1, nsChan, :) = permute(NScoefMat(nsIndOrig, :), [3 1 2]);
coefMat(2, wfsChan, :) = permute(WFScoefMat*corrVol(wfsIndOrig, :), [3 1 2]);
coefMat(3, nsChan, :) = permute(NScoefMat(nsIndOrig, :), [3 1 2]);
coefMat(3, wfsChan, :) = permute(WFScoefMat*corrVol(wfsIndOrig, :), [3 1 2]);

pulseLimits = [1 3; 4 6; 7 9];

sampleRate = 44100/4;
signalFunction = @(startSample, endSample) pulseCoefMat2signal(coefMat, pulseLimits,...
    freqs, sampleRate, startSample, endSample, 'type_pulseLimits', 'time');

obj.WFSToolObj.reproduceSignalFunction(signalFunction, sampleRate);
recFieldCorrVol = obj.WFSToolObj.recorded;

% ---- Optimized volume and noise source location ----
% No completado

%% Comparison of measures and estimations. They should be similar.

% Get the received signal frequency spectrum

% Compare it with the theoretical one


%% View of results


