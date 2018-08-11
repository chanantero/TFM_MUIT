%% Experiment 16.

% Created on 10/08/2018.

% Measures on the GTAC listening room. Reproduction and recording.

%% Preamble
pathSetUp;

imagesPath = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\';

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
    coefMat(inds) = exp(1i*rand*2*pi);
end

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

% Load frequency response test signal parameters
IDload = '2018-08-11_19-02-13';
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

% Reproduce
filename = [dataPathName, 'acousticPathReprodSignal_', IDload, '.wav'];
repRecObj = reproductorRecorder;
repRecObj.setProps('mode', originType('file'), 1);
repRecObj.setProps('enableProc', false);
% filename = 'C:\Users\Rubén\Music\Salsa\Flor Pálida - Marc Anthony.mp3';
repRecObj.setProps('audioFileName', filename, 1);
repRecObj.executeOrder('play');

% Reproduce in simulation
[y, Fs] = audioread(filename, 'native');

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
recPositions = [centerX, centerY, 0];

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
NSpositions = [0 0 0];

SetupParametersScript
AcousticPathCalculationScript

obj.domain = 'time';

WFSacPathIR = permute(WFS_IR(:, :, :, 1), [1, 3, 2]);
obj.setAcousticPaths('WFS', WFSacPathIR);

fragLength = 2^12; % fragment length in samples
numSamp = size(y, 1);
numFrag = ceil(numSamp/fragLength); % number of fragments
previousSufix = zeros(numMicro, numSampIR - 1);
recSignal = zeros(numMicro, numSamp); % Recorded signal. (numMicros x numSamples)
for f = 1:numFrag
    fprintf('%d/%d\n', f, numFrag);
    selSamp = 1 + fragLength*(f - 1):min(fragLength*f, numSamp); % The min() is for the last iteration
    frag = [y(selSamp, :)', zeros(numChannels, numSampIR-1)];
    obj.WFSToolObj.WFSarrayCoefficient = frag;
    obj.WFSToolObj.simulTheo.updateField();
    recFrag = obj.WFSToolObj.simulField;
    recFrag(:, 1:numSampIR-1) = recFrag(:, 1:numSampIR-1) + previousSufix; % Add previous sufix
    previousSufix = recFrag(:, end - (numSampIR - 1) + 1:end); % Set sufix for next loop iteration
    recSignal(:, selSamp) = recFrag(:, 1:length(selSamp)); % Save the signal without the sufix. We use length(selSamp) and not fragLength for the last iteration
end

% Crop signal with an adequate duration
numMicros = size(recSignal, 1);

margin = 1;

numPulses = size(pulseLimits, 1);
FR = zeros(numMicros, numChannels, numFreqs);
for p = 1:numPulses
    start = pulseLimits(p, 1) + margin;
    ending = pulseLimits(p, 2) - margin;
    
    startIndex = floor(start*sampleRate) + 1;
    endIndex = floor(ending*sampleRate) + 1;
    
    for m = 1:numMicros
        fragment = recSignal(m, startIndex:endIndex);
        X = freqz(fragment, 1, freqs, sampleRate); % (1 x numFreqs)
        
        % Cómo guardar la información si de cada canal tenemos solo
        % información de algunas frecuencias?
        % Idea, da igua, llena la matriz y los valores a 0 se entiende que
        % han de interpolarse
        
        coefMatFlag = permute(coefMat(p, :, :) ~= 0, [3 2 1]);
        for c = 1:numChannels
            FR(m, c, coefMatFlag(:, c)) = X(coefMatFlag(:, c));
        end       
    end
end

% Check the result is equal to the old form of doing it. All must be
% simulated (by now) of course.
