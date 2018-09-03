%% Experiment 16.

% Created on 10/08/2018.

% Measures on the GTAC listening room. Reproduction and recording.

%% Preamble
globalPath = './';
paths = genpath(globalPath);
addpath(paths);

dataPathName = [globalPath, 'Data\'];
ID = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

%% System set up.
obj = SimulationController;
obj.ax.Parent.HandleVisibility = 'off';
obj.ax.HandleVisibility = 'off';
obj.WFSToolObj.fig.HandleVisibility = 'off';

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

% ---- Generate signal and save it into a file ----
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

% Next line produces underruns
% obj.WFSToolObj.reproduceSignalFunction(signalFunction, sampleRate);
% So, it is better to read from file

% Write into file
audWriterObj = dsp.AudioFileWriter;
audWriterObj.Filename = [dataPathName, 'test', ID, '.wav'];
audWriterObj.SampleRate = sampleRate;
audWriterObj.DataType = 'single';

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

% ---- Reproduce ----
IDload = 'Lab_23-08-2018';
filename = [dataPathName, 'test', IDload, '.wav'];
info = audioinfo(filename);
sampleRate = info.SampleRate;

samplesPerFrame = sampleRate*2;
numRecChann = 1;
numFrames = ceil(info.TotalSamples/samplesPerFrame);

signProv = signalProvider;
signProv.FileName = filename;
signProv.SamplesPerFrame = samplesPerFrame; % The maximum an audioDeviceWriter allows

% Reproduce with audioDeviceWriter object
deviceWriter = audioDeviceWriter;
deviceWriter.SampleRate = sampleRate;
deviceWriter.Driver = 'ASIO';
deviceWriter.Device = 'MOTU PCI ASIO';
% deviceWriter.ChannelMappingSource = 'Property';

count = 0;
while ~isDone(signProv)
    x = step(signProv);
    step(deviceWriter, x);
    count = count + 1;
    fprintf('%d\n', count);
end
release(signProv)
release(deviceWriter)

%% Test that the OFDM and Chirp acoustic path estimation are equivalent

% ---- OFDM Channel Estimation ----
% Generate signal and write it into a file
numTotalChan = 2;
numSimultChan = 2; % Number of simultaneous channels reproducing
freqRange = [20, 1200]; % Hz. It's a user preference, but it's not going to be exactly respectd
freqStep = 1; % Hz
pulseDurationInPeriods = 60;
silenceBetweenPulses = 3;
sampleRate = 22050;

filename = [dataPathName, 'signalStructTestOFDM'];
OFDMchannelEstimationSignal(numTotalChan, numSimultChan, freqRange, freqStep, pulseDurationInPeriods, silenceBetweenPulses, sampleRate, 'filename', filename);

% Load frequency response test signal parameters
IDload = 'signalStructTestOFDM';
signalStructLoaded = load([dataPathName, IDload, '.mat']);
pulseLimits = signalStructLoaded.pulseLimits;
freqs = signalStructLoaded.freqs;
filename = [dataPathName, IDload, '.wav'];
info = audioinfo(filename);
sampleRate = info.SampleRate;
coefMat = signalStructLoaded.coefMat;

% Reproduce and record
samplesPerFrame = sampleRate*2;
playRecExt = audioPlayerRecorderExtended;
playRecExt.mode = originType('file');
playRecExt.filename = filename;
playRecExt.SampleRate = sampleRate;
playRecExt.SamplesPerFrame = samplesPerFrame;
location = 'laptop'; % laboratory/laptop
switch location
    case 'laboratory'
        playRec.Device = 'MOTU PCI ASIO';
    case 'laptop'
        playRec.Device = 'Default';
        numRecChann = 1;
end
playRecExt.RecorderChannelMapping = 1:numRecChann;
playRec.PlayerChannelMapping = [1 2];
playRecExt.playAndRecord();
y = playRecExt.recSignal;  

% ax = axes(figure, 'NextPlot', 'Add');
% plot(ax, x)
% plot(ax, y)
% sound(y(44100*3+1:end), sampleRate)

% Correct offset between the transmitted and received signals
numSamp = size(y, 1);
t = (0:numSamp-1)/sampleRate;
ind = find(t < pulseLimits(1,2)); % Consider only the first pulse
x = audioread(filename, [1, ind(end)], 'native');
cor = xcorr(mean(x, 2), y(ind, 1));
ax = axes(figure);
plot(ax, cor)
centerCor = (length(cor) + 1)/2;
[~, indMax] = max(cor);
numShiftFrames = round((indMax - centerCor)/samplesPerFrame);
recSignal = shiftArray(y, numShiftFrames*samplesPerFrame);

% Estimate channel
Hofdm = OFDMchannelEstimation( filename, recSignal, coefMat, pulseLimits, freqs, sampleRate);

% Interpolate
FRint = zeros(size(Hofdm));
numMicros = size(recSignal, 2);
numChannels = size(x, 2);
for m = 1:numMicros
    for c = 1:numChannels
        nonZero = Hofdm(m, c, :) ~= 0; % recSpec(m, c, :) ~= 0;
        magn = interp1(freqs(nonZero), abs(squeeze(Hofdm(m, c, nonZero))), freqs, 'linear');
        phase = interp1(freqs(nonZero), unwrap(angle(squeeze(Hofdm(m, c, nonZero)))), freqs, 'linear');
        FRint(m, c, :) = magn.*exp(1i*phase);
    end
end

% ax = axes(figure);
% plot(ax, freqs, abs(squeeze(FRint)))


% ---- Chirp estimation ----

% Generate chirp signal
durSign = 30; % Duration of tone for time processing
t = (0:ceil(durSign*sampleRate)-1)/sampleRate;
chirpSignal = 0.15*chirp(t, freqRange(1), durSign, freqRange(2));
prefixDuration = 2;
prefixNumSamples = ceil(prefixDuration*sampleRate);
chirpSignal = [zeros(1, prefixNumSamples), chirpSignal, zeros(1, prefixNumSamples )]';
t = (0:length(chirpSignal)-1)/sampleRate;

samplesPerFrame = sampleRate*2;
playRecExt = audioPlayerRecorderExtended;
playRecExt.PlayerChannelMapping = 1;
playRecExt.mode = originType('custom');
playRecExt.customSignal = chirpSignal;
playRecExt.SampleRate = sampleRate;
playRecExt.SamplesPerFrame = samplesPerFrame;
playRecExt.Device = 'Default';
numRecChann = 1;
playRecExt.RecorderChannelMapping = 1:numRecChann;
playRecExt.playAndRecord();
recSignal1 = playRecExt.recSignal;

playRecExt.PlayerChannelMapping = 2;
playRecExt.playAndRecord();
recSignal2 = playRecExt.recSignal;

specTransChirp = freqz(chirpSignal, 1, freqs, sampleRate);
specRecChirp1 = freqz(recSignal1, 1, freqs, sampleRate);
specRecChirp2 = freqz(recSignal2, 1, freqs, sampleRate);

% ax = axes(figure, 'NextPlot', 'Add');
% plot(ax, freqs, abs(specRecChirp1))
% plot(ax, freqs, abs(specTransChirp))

Hchirp1 = specRecChirp1./specTransChirp;
Hchirp2 = specRecChirp2./specTransChirp;
Hchirp = [Hchirp1; Hchirp2];

ax = axes(figure, 'NextPlot', 'Add');
lchirp = plot(ax, freqs, abs(Hchirp));
lofdm = plot(ax, freqs, abs(squeeze(FRint)));

lchirp(1).LineStyle = ':';
lchirp(2).LineStyle = ':';

% Da lo mismo aproximadamente!! :D eso es que el m�todo funciona

% Different amplitude
dur = 10;
t = (0:1/sampleRate:dur)';
amp = linspace(0, 1, length(t))';
freq = 440;
signal = amp.*cos(2*pi*freq*t);

playRecExt.PlayerChannelMapping = 1;
playRecExt.customSignal = signal;
playRecExt.playAndRecord();
recSignal = playRecExt.recSignal;

ax = axes(figure, 'NextPlot', 'Add');
plot(ax, signal)
plot(ax, 15*recSignal)

%% Generate frequency response signal

% In order to measure the acoustic paths, we are going to use the method of
% orthogonal frequecies.
% We are going to reproduce evenly spaced tones in a given frequency range.
% However, every loudspeaker is going to reproduce different frequencies,
% so, in the reception of signals, we can, by spectral analysis,
% differenciate the signal from the different loudspeakers.

% Generate signal and write it into a file
numTotalChan = 96;
numSimultChan = 4; % Number of simultaneous channels reproducing
freqRange = [20, 1200]; % Hz. It's a user preference, but it's not going to be exactly respectd
freqStep = 1; % Hz
pulseDurationInPeriods = 20;
silenceBetweenPulses = 3;
sampleRate = 44100;

filename = [dataPathName, 'acousticPathReprodSignal_', ID, '.wav'];
OFDMchannelEstimationSignal(numTotalChan, numSimultChan, freqRange, freqStep, pulseDurationInPeriods, silenceBetweenPulses, sampleRate, 'filename', filename);

%% Reproduce and record frequency response test signal

% Load frequency response test signal parameters
IDload = 'Lab_29-08-2018';
signalStructLoaded = load([dataPathName, 'signalStruct_', IDload, '.mat']);
numChannels = 96;
pulseLimits = signalStructLoaded.pulseLimits;
freqs = signalStructLoaded.freqs;
numFreqs = length(freqs);
coefMat = signalStructLoaded.coefMat;
filename = [dataPathName, 'acousticPathReprodSignal_', IDload, '.wav'];
info = audioinfo(filename);
sampleRate = info.SampleRate;

% Reproduce and record
simulatedReproduction = true;
if ~simulatedReproduction
    audioPlayerRecorderEnabled = strcmp(version('-release'), '2018a');
    if ~audioPlayerRecorderEnabled
        % Reproduce
        repRecObj = reproductorRecorder;
        repRecObj.setProps('mode', originType('file'), 1);
        repRecObj.setProps('enableProc', false);
        % filename = 'C:\Users\Rub�n\Music\Salsa\Flor P�lida - Marc Anthony.mp3';
        repRecObj.setProps('audioFileName', filename, 1);
        repRecObj.executeOrder('play');
        recSignal = repRecObj.recorded{1}.';
    else
        % In case it doesn't work: DebuggingGTAC A).
        samplesPerFrame = sampleRate*2;
        playRecExt = audioPlayerRecorderExtended;
        playRecExt.mode = originType('file');
        playRecExt.filename = filename;
        playRecExt.SampleRate = sampleRate;
        playRecExt.SamplesPerFrame = samplesPerFrame;
        location = 'laptop'; % laboratory/laptop
        switch location
            case 'laboratory'
                playRec.Device = 'MOTU PCI ASIO';
            case 'laptop'
                playRec.Device = 'Default';
                numRecChann = 1;
        end
        playRecExt.RecorderChannelMapping = 1:numRecChann;
        playRec.PlayerChannelMapping = [1 2];
        playRecExt.playAndRecord();
        y = playRecExt.recSignal;
            
        % Adjust possible delay between tranmitted and recorded signal. It
        % is typical that the delay is of one frame.
        numSamp = size(y, 1);
        t = (0:numSamp-1)/sampleRate;
        ind = t < pulseLimits(2,1); % Consider only the first pulse
        cor = xcorr(mean(x(ind, :), 2), y(ind, 1));
        ax = axes(figure);
        plot(ax, cor)
        centerCor = (length(cor) + 1)/2;
        [~, indMax] = max(cor);
        numShiftFrames = round((indMax - centerCor)/samplesPerFrame); 
        recSignal = shiftArray(y, numShiftFrames*samplesPerFrame);
        recSignal = recSignal';
        
%         ax = axes(figure, 'NextPlot', 'Add');
%         plot(ax, x(1:samplesPerFrame*5, 1))
%         plot(ax, y(1:samplesPerFrame*5, 1))
%         plot(ax, recSignal(1, 1:samplesPerFrame*5))
    end
else
    % Reproduce in simulation    
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
    info = audioinfo(filename);
    numSamp = info.TotalSamples;
    numFrag = ceil(numSamp/fragLength); % number of fragments
    previousSufix = zeros(numMicro, numSampIR - 1);
    recSignal = zeros(numMicro, numSamp); % Recorded signal. (numMicros x numSamples)
    for f = 1:numFrag
        fprintf('%d/%d\n', f, numFrag);
        firstSamp = 1 + fragLength*(f - 1);
        lastSamp = min(fragLength*f, numSamp); % The min() is for the last iteration
        fragSignal = audioread(filename, [firstSamp, lastSamp]);
        frag = [fragSignal', zeros(numChannels, numSampIR-1)];
        obj.WFSToolObj.WFSarrayCoefficient = zeros(obj.numWFS, size(frag, 2));
        obj.WFSToolObj.WFSarrayCoefficient(obj.WFSToolObj.loudspeakerMapping(1).originInd, :) = frag(wfsChan, :);
        obj.WFSToolObj.noiseSourceCoefficient_complete = zeros(2, size(frag, 2));
        obj.WFSToolObj.noiseSourceCoefficient_complete(obj.WFSToolObj.loudspeakerMapping(2).originInd, :) = frag(nsChan, :);
        obj.WFSToolObj.simulTheo.updateField();
        recFrag = obj.WFSToolObj.simulField;
        recFrag(:, 1:numSampIR-1) = recFrag(:, 1:numSampIR-1) + previousSufix; % Add previous sufix
        previousSufix = recFrag(:, end - (numSampIR - 1) + 1:end); % Set sufix for next loop iteration
        selSamp = firstSamp:lastSamp;
        recSignal(:, selSamp) = recFrag(:, 1:length(selSamp)); % Save the signal without the sufix. We use length(selSamp) and not fragLength for the last iteration
    end
end

FR = OFDMchannelEstimation( filename, recSignal', coefMat, pulseLimits, freqs, sampleRate);

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
%     fig = figure;
%     ax = axes(fig);
%     histogram(ax, abs(rel(recFlag)))
%     histogram(angle(rel(recFlag)))
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
FRint = zeros(size(FR));
for m = 1:numMicros
    for c = 1:numChannels
        nonZero = FR(m, c, :) ~= 0; % recSpec(m, c, :) ~= 0;
        magn = interp1(freqs(nonZero), abs(squeeze(FR(m, c, nonZero))), freqs, 'linear');
        phase = interp1(freqs(nonZero), unwrap(angle(squeeze(FR(m, c, nonZero)))), freqs, 'linear');
        FRint(m, c, :) = magn.*exp(1i*phase);
    end
end

% save([dataPathName, 'acousticPathRecordSignal_', ID, '.mat'], 'IDload', 'recSignal', 'FR')
%% Estimation of results based on measured acoustic path responses

% ---- No optimization ----
% Set parameters
% Calculate noise source position based on distances measured in the lab
w = 0.2;
d = 0.184;
L = 1.034;
alpha = 45;
x = (1 + cosd(alpha))*8*0.18 + w*sind(alpha) + d*cosd(alpha) + L*sind(alpha);
y = 0 - w*cosd(alpha) + d*sind(alpha) - L*cosd(alpha);
NSpositions = [x y 0]; % Assumed real position

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
    
    sSimul = s';

    [sExtTheo, corrFactInd, corrFactGlob, gainInd, gainGlob, corrFactAver, gainAver] =...
        SimulationController.addCancellationParametersToStructure(sSimul);
    
    ax = axes(figure);
    plot(ax, freqs, 10*log10(gainAver));
    plot(ax, freqs, abs(corrFactAver))
    plot(ax, freqs, unwrap(angle(corrFactAver)))
    
else

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
corrFactComp = -fieldNSaux./fieldWFSaux;
plot(freqs(sel), abs(corrFactComp))

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
sampleRate = 44100;
t = (0:ceil(durSign*sampleRate)-1)/sampleRate;
NSsignal = chirp(t, min(freqs), durSign, max(freqs) + 50);
prefixDuration = 2;
prefixNumSamples = ceil(prefixDuration*sampleRate);
NSsignal = [zeros(1, prefixNumSamples), NSsignal, zeros(1, prefixNumSamples )];
t = (0:length(NSsignal)-1)/sampleRate;

% Frequency filters
magnFiltOrder = 2^12;
hilbertFiltOrder = 2^12;
[freqFilter, freqFiltDelay] = getFrequencyFilter(magnFiltOrder, hilbertFiltOrder, sampleRate, 'analytical', true);
obj.WFSToolObj.freqFilter = freqFilter;
% Debug
% freqFilterResp = freqz(freqFilter, 1, freqs, sampleRate);
% ax = axes(figure);
% plot(ax, freqs, abs(freqFilterResp), freqs, sqrt(freqs/340))
% plot(ax, freqs, rad2deg(wrapToPi(angle(freqFilterResp) + freqFiltDelay/sampleRate*2*pi*freqs)))

% Calculate WFS signals
obj.domain = 'time';
obj.NScoef = NSsignal;
obj.NSVcoef = -NSsignal;

frequencyCorrection = true;
obj.NSposition = NSpositions;
obj.WFSToolObj.frequencyCorrection = true;
obj.WFSToolObj.attenuationType = 'Ruben'; attenuationType = 'Ruben';
obj.WFSToolObj.updateFiltersWFS();
obj.WFSToolObj.automaticLengthModification = false;
obj.WFSToolObj.WFScalculation();

% Reproduce and record:
if ~simulatedReproduction
    numSamp = length(NSsignal);
    customSignal = zeros(numSamp, numChannels);
    wfsIndOrig = obj.WFSToolObj.loudspeakerMapping(1).originInd;
    wfsIndDest = obj.WFSToolObj.loudspeakerMapping(1).destinationInd;
    nsIndDest = obj.WFSToolObj.loudspeakerMapping(2).destinationInd;
    customSignal(:, nsIndDest) = NSsignal(:);
    customSignal(:, wfsIndDest) = obj.WFScoef(wfsIndOrig, :).';
    
    % In case audioPlayerRecorder doesn't work, consult debuggingGTAC B)
%     repRecObj = reproductorRecorder;
%     repRecObj.setProps('mode', originType('custom'), 1);
%     repRecObj.setProps('enableProc', false);
%     repRecObj.setProps('Fs_player', sampleRate, 1);
%     repRecObj.setProps('Fs_recorder', sampleRate, 1);
    
    samplesPerFrame = sampleRate*2;
    numFrames = ceil(numSamp/samplesPerFrame);

    playRec = audioPlayerRecorder;
    playRec.SampleRate = sampleRate;
    playRec.Device = 'MOTU PCI ASIO';
    numRecChann = 2;
    playRec.RecorderChannelMapping = 1:numRecChann;
    
    signProv = signalProvider;
    signProv.mode = originType('custom');
    signProv.SamplesPerFrame = samplesPerFrame;
    
% ---- No optimization ----
    
    % - Only noise source
    aux = customSignal;
    aux(:, wfsIndDest) = 0;
%     repRecObj.setProps('customSignal', aux, 1);
%     repRecObj.executeOrder('play');
%     recSignalNS = repRecObj.recorded{1}.';
    signProv.customSignal = aux;
    recSignalNS = zeros(numSamp, numRecChann);
    ind = 1-samplesPerFrame:0;
    for fr = 1:numFrames
        x = step(signProv);
        ind = ind + samplesPerFrame;
        recSignalNS(ind, :) = step(playRec, x);
    end
    release(signProv)
    release(playRec)
    
    cor = xcorr(mean(aux, 2), recSignalNS(:, 1));
    ax = axes(figure);
    plot(ax, cor)
    plot(ax, aux(:, 1))
    hold on
    plot(ax, recSignalNS(:, 1))
    centerCor = (length(cor) + 1)/2;
    [~, indMax] = max(cor);
    numShiftFrames = round((indMax - centerCor)/samplesPerFrame); 
    recSignalNS = shiftArray(recSignalNS, numShiftFrames*samplesPerFrame);
    recSignalNS = recSignalNS';
    
    % - Only WFS
    aux = customSignal;
    aux(:, nsIndDest) = 0;
%     repRecObj.setProps('customSignal', aux);
%     repRecObj.executeOrder('play');
%     recSignalWFS = repRecObj.recorded{1}.';
    signProv.customSignal = aux;
    recSignalWFS = zeros(numSamp, numRecChann);
    ind = 1-samplesPerFrame:0;
    for fr = 1:numFrames
        ind = ind + samplesPerFrame;
        x = step(signProv);
        recSignalWFS(ind, :) = step(playRec, x);
    end
    release(signProv)
    release(playRec)
    
    cor = xcorr(mean(aux, 2), mean(recSignalWFS, 2));
    centerCor = (length(cor) + 1)/2;
    [~, indMax] = max(cor);
    numShiftFrames = round((indMax - centerCor)/samplesPerFrame); 
    recSignalWFS = shiftArray(recSignalWFS, numShiftFrames*samplesPerFrame);
    recSignalWFS = recSignalWFS';  

    % - All
%     repRecObj.setProps('customSignal', customSignal);
%     repRecObj.executeOrder('play');
%     recSignal = repRecObj.recorded{1}.';
    signProv.customSignal = customSignal;
    recSignal = zeros(numSamp, numRecChann);
    ind = 1-samplesPerFrame:0;
    for fr = 1:numFrames
        ind = ind + samplesPerFrame;
        x = step(signProv);
        recSignal(ind, :) = step(playRec, x);
    end
    release(signProv)
    release(playRec)
    
    cor = xcorr(mean(aux, 2), mean(recSignal, 2));
    centerCor = (length(cor) + 1)/2;
    [~, indMax] = max(cor);
    numShiftFrames = round((indMax - centerCor)/samplesPerFrame); 
    recSignal = shiftArray(recSignal, numShiftFrames*samplesPerFrame);
    recSignal = recSignal';  
    
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
    c = 340;
        
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
%     ax = axes(figure);
%     plot(ax, t, NSsignal, t, obj.WFScoef(89, :))
%     nsSignalFreq = freqz(NSsignal, 1, freqs, sampleRate);
%     wfsSignalFreq = freqz(obj.WFScoef(89, :), 1, freqs, sampleRate);
%     ax = axes(figure);
%     dist = calcDistances(obj.NSRposition, obj.WFSposition(89,:));
%     cosAlpha = dot((obj.WFSposition(89,:) - obj.NSRposition), [0 1 0]);
%     plot(ax, freqs, abs(wfsSignalFreq./nsSignalFreq))
%     plot(ax, freqs, rad2deg(wrapToPi(angle(wfsSignalFreq./nsSignalFreq) + dist/340*2*pi*freqs)))
%     plot(ax, freqs, angle(wfsSignalFreq./nsSignalFreq) )
%     
%     % size(obj.filtersWFS_IR) (numWFS x numNS x numSamp)
%     filt = squeeze(obj.WFSToolObj.filtersWFS_IR(89, 2, :));
%     plot(ax, filt)
%     filtSpec = freqz(filt, 1, freqs, sampleRate);
%     plot(ax, freqs, angle(filtSpec))
%     
%     wfsSignalTheo = fftfilt(filt, NSsignal);
%     plot(ax, t, NSsignal, t, wfsSignalTheo);
    
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

% save([dataPathName, 'recSignals_', ID, '.mat'], 'sampleRate', 'NSsignal', 'recSignalNS', 'recSignalWFS', 'recSignal')%, ...
    %'recSignalNScorrVol', 'recSignalWFScorrVol', 'recSignalCorrVol');

%% Comparison of measures and estimations. They should be similar.
% IDload = 'Lab_29-08-2018_corregido';
% load([dataPathName, 'recSignals_', ID, '.mat'])

numMicro = size(recSignal, 1);
corrVolTime = zeros(numMicro, 1);
for m = 1:numMicro
    corrVolTime(m) = -recSignalNS(m,:)'\recSignalWFS(m,:)';
end

numSamp = size(recSignal, 2);
recSignalNScorrVol = recSignalNS;
recSignalWFScorrVol = recSignalWFS.*repmat(corrVolTime, [1, numSamp]);
recSignalCorrVol = recSignalNScorrVol + recSignalWFScorrVol;

% Get the received signal frequency spectrum
oper = @(x) freqz(x, 1, freqs, sampleRate);

NSspec = freqz(NSsignal, 1, freqs, sampleRate);
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
% [sExtSimulCorrVol, corrFactIndSimulCorrVol, corrFactGlobSimulCorrVol, gainIndSimulCorrVol, gainGlobSimulCorrVol, corrFactAverSimulCorrVol, gainAverSimulCorrVol] =...
%         SimulationController.addCancellationParametersToStructure(sSimulCorrVol);

ax = axes(figure, 'NextPlot', 'Add');
plot(ax, freqs, 10*log10(gainAver))
plot(ax, freqs, 10*log10(gainAverSimul))
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = 'Average gain (dB)';

ax = axes(figure, 'NextPlot', 'Add');
plot(ax, freqs, abs(corrFactIndSimul(:, 1, 1)))
plot(ax, freqs, abs(corrFactInd(:, 1, 1)))
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = '|\Psi|';

% ax = axes(figure, 'NextPlot', 'Add');
% plot(ax, freqs, 10*log10(gainAverCorrVol))
% plot(ax, freqs, 10*log10(gainAverSimulCorrVol))
% ax.XLabel.String = 'Frequency (Hz)';
% ax.YLabel.String = 'Average gain (dB)';

% % Debug
% % The received signal spectrum from the noise source is the same
% NSspec = freqz(NSsignal, 1, freqs, sampleRate);
% resp = recNS./repmat(NSspec, 2, 1);
% recNSsimul = [sSimul.recNScoef];
% ax = axes(figure);
% plot(ax, freqs, abs(resp), freqs, abs(recNSsimul))
% 
% % The received signal spectrum from the WFS array is the same? No, but it's
% % similar. I think the differece is due to the fact that FRint is not
% % exactly as the simulated acoustic paths, since it is an interpolation. Or
% % maybe it is because the frequency response of the frequency filter is not
% % ideal
% resp = recWFS./repmat(NSspec, 2, 1);
% recWFSsimul = [sSimul.recWFScoef];
% ax = axes(figure);
% plot(ax, freqs, abs(resp), freqs, abs(recWFSsimul)) % plot(ax, freqs, angle(resp), freqs, angle(recNSsimul))
% 
% ax = axes(figure);
% plot(ax, t, recSignalNS(1,:), t, recSignalWFS(1,:), t, recSignal(1,:))

%% View of results
numSamp = size(NSsignal, 2);
tNS = (0:numSamp - 1)/sampleRate;
numSampRec = size(recSignal, 2);
tRec = (0:numSampRec - 1)/sampleRate;

% The transmitted signal by the noise source is:
    % Time representation
    ax = axes(figure);
    plot(ax, tNS, NSsignal)
    ax.XLabel.String = 'Time (s)';
    ax.YLabel.String = '$\signal[ns][time]$'; ax.YLabel.Interpreter = 'latex';
%     options.TickLabels2Latex = false;
%     Plot2LaTeX(ax.Parent, [imagesPath, 'Experiment16_NSsignalTime'], options)

    % Frequency representation
    ax = axes(figure);
    plot(ax, freqs, abs(NSspec));
    ax.XLabel.String = 'Frequency (Hz)';
    ax.YLabel.String = '$\signal[ns][frequency]$'; ax.YLabel.Interpreter = 'latex';
%     Plot2LaTeX(ax.Parent, [imagesPath, 'Experiment16_NSsignalFreq'], options)
    
% The received signals from the noise source in both microphones are:
    % Time representation
    ax = axes(figure);
    plot(ax, tRec, recSignalNS.')
    ax.XLabel.String = 'Time (s)';
    ax.YLabel.String = '$\Field[ns][time]$'; ax.YLabel.Interpreter = 'latex';
    Plot2LaTeX(ax.Parent, [imagesPath, 'Experiment16_recNSTime'], options)

    % Frequency representation
    ax = axes(figure);
    plot(ax, freqs, abs(recNS).')
    ax.XLabel.String = 'Frequency (Hz)';
    ax.YLabel.String = '$\Field[ns][frequency]$'; ax.YLabel.Interpreter = 'latex';
%     Plot2LaTeX(ax.Parent, [imagesPath, 'Experiment16_recNSFreq'], options)

% The received total signals are:
    % ---- No optimization ----
        % Time representation
        ax = axes(figure);
        plot(ax, tRec, recSignal.')
        ax.XLabel.String = 'Time (s)';
        ax.YLabel.String = '$\Field[total][time]$'; ax.YLabel.Interpreter = 'latex';
%         Plot2LaTeX(ax.Parent, [imagesPath, 'Experiment16_recTime'], options)

        % Frequency representation
        ax = axes(figure);
        plot(ax, freqs, abs(rec).')
        ax.XLabel.String = 'Frequency (Hz)';
        ax.YLabel.String = '$\Field[total][frequency]$'; ax.YLabel.Interpreter = 'latex';
%         Plot2LaTeX(ax.Parent, [imagesPath, 'Experiment16_recFreq'], options)

    % ---- Volume optimization ----    
        % Time representation
        ax = axes(figure);
        plot(ax, tRec, recSignalCorrVol.')
        ax.XLabel.String = 'Time (s)';
        ax.YLabel.String = '$\Field[total][time]$'; ax.YLabel.Interpreter = 'latex';
%         Plot2LaTeX(ax.Parent, [imagesPath, 'Experiment16_recCorrVolTime'], options)

        % Frequency representation
        ax = axes(figure);
        plot(ax, freqs, abs(recCorrVol).')
        ax.XLabel.String = 'Frequency (Hz)';
        ax.YLabel.String = '$\Field[total][frequency]$'; ax.YLabel.Interpreter = 'latex';
%         Plot2LaTeX(ax.Parent, [imagesPath, 'Experiment16_recCorrVolFreq'], options)

% The received total signals in comparison with the receved noise source signals
    % ---- No optimization ----
        % Time representation
        ax = axes(figure);
        plot(ax, tRec, recSignalNS(1,:).', tRec, recSignal(1,:).')
        ax.XLabel.String = 'Time (s)';
        legend(ax, {'$\Field[ns][time]$', '$\Field[total][time]$'})
%         Plot2LaTeX(ax.Parent, [imagesPath, 'Experiment16_recAndrecNStime_1'], options)

        ax = axes(figure);
        plot(ax, tRec, recSignalNS(2,:).', tRec, recSignal(2,:).')
        ax.XLabel.String = 'Time (s)';
        legend(ax, {'$\Field[ns][time]$', '$\Field[total][time]$'})
%         Plot2LaTeX(ax.Parent, [imagesPath, 'Experiment16_recAndrecNStime_2'], options)

        % Frequency representation
        ax = axes(figure);
        plot(ax, freqs, abs(recNS(1,:)).', freqs, abs(rec(1,:)).')
        ax.XLabel.String = 'Frequency (Hz)';
        legend(ax, {'$\Field[ns][frequency]$', '$\Field[total][frequency]$'})
%         Plot2LaTeX(ax.Parent, [imagesPath, 'Experiment16_recAndrecNSfreq_1'], options)

        ax = axes(figure);
        plot(ax, freqs, abs(recNS(2,:)).', freqs, abs(rec(2,:)).')
        ax.XLabel.String = 'Frequency (Hz)';
        legend(ax, {'$\Field[ns][frequency]$', '$\Field[total][frequency]$'})
%         Plot2LaTeX(ax.Parent, [imagesPath, 'Experiment16_recAndrecNSfreq_2'], options)

% ---- Volume optimization ----
        % Time representation
        ax = axes(figure);
        plot(ax, tRec, recSignalNScorrVol(1,:).', tRec, recSignalCorrVol(1,:).')
        ax.XLabel.String = 'Time (s)';
        legend(ax, {'$\Field[ns][time]$', '$\Field[total][time]$'})
%         Plot2LaTeX(ax.Parent, [imagesPath, 'Experiment16_recAndrecNStime_1'], options)

        ax = axes(figure);
        plot(ax, tRec, recSignalNScorrVol(2,:).', tRec, recSignalCorrVol(2,:).')
        ax.XLabel.String = 'Time (s)';
        legend(ax, {'$\Field[ns][time]$', '$\Field[total][time]$'})
%         Plot2LaTeX(ax.Parent, [imagesPath, 'Experiment16_recAndrecNStime_2'], options)

        % Frequency representation
        ax = axes(figure);
        plot(ax, freqs, abs(recNScorrVol(1,:)).', freqs, abs(recCorrVol(1,:)).')
        ax.XLabel.String = 'Frequency (Hz)';
        legend(ax, {'$\Field[ns][frequency]$', '$\Field[total][frequency]$'})
%         Plot2LaTeX(ax.Parent, [imagesPath, 'Experiment16_recAndrecNSfreq_1'], options)

        ax = axes(figure);
        plot(ax, freqs, abs(recNScorrVol(2,:)).', freqs, abs(recCorrVol(2,:)).')
        ax.XLabel.String = 'Frequency (Hz)';
        legend(ax, {'$\Field[ns][frequency]$', '$\Field[total][frequency]$'})
%         Plot2LaTeX(ax.Parent, [imagesPath, 'Experiment16_recAndrecNSfreq_2'], options)


% Average gain for no optimization and for optimization of volume
ax = axes(figure, 'NextPlot', 'Add');
plot(ax, freqs, 10*log10(gainAver), freqs, 10*log10(gainAverCorrVol))
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = 'Average gain (dB)';
legend(ax, 'No optimization', 'Volume correction')
% printfig(ax.Parent, imagesPath, 'Experiment16_averGainComparison', 'eps')