%% Reproduction tests in laptop

%% Preamble
pathSetUp;

imagesPath = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\';

dataPathName = [globalPath, 'Data\'];
ID = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

%% System set up.
obj = SimulationController;

obj.NSposition = [3.35 -0.2 0]; % Assumed real position
obj.amplitude = 1;
obj.phase = 0;
obj.frequency = 440;

% Testing for laptop
obj.WFSToolObj.setNumWFSarraySources(2);
obj.WFSToolObj.setNumReceivers(1);

% Set noise source channel through the GUI or use next line
NSchan = 1;
obj.WFSToolObj.noiseSourceChannelMapping(1) = NSchan;

%% Test that loudspeakers reproduce in an adequate way
numWFS = 96;
coefMat = zeros(numWFS, 2, 8); %zeros(numWFS, numWFS, 8);
semiToneShift = [-9 -7 -5 -4 -2 0 2 3];
freqLA = 440;
freqs = freqLA .* 2.^(semiToneShift/12);
prog = repmat([(1:7), (8:-1:2)], 1, ceil(numWFS/14));
indFreq = prog(1:numWFS);
indChan = mod(1:numWFS, 2)+1;%1:numWFS;
indPulse = 1:numWFS;
inds = sub2ind([numWFS, 2, 8], indPulse, indChan, indFreq);
coefMat(inds) = 1;

pulseDur = 1.5; % seconds
pulseSolap = 0.5; % seconds
pulseStart = (0:numWFS-1)'*(pulseDur - pulseSolap);
pulseEnd = pulseStart + pulseDur;
pulseLimits = [pulseStart, pulseEnd];

sampleRate = 44100;
signalFunction = @(startSample, endSample) pulseCoefMat2signal(coefMat, pulseLimits,...
    freqs, sampleRate, startSample, endSample, 'type_pulseLimits', 'time');

obj.WFSToolObj.reproduceSignalFunction(signalFunction, 44100);

%% Test for theoric frequency response
obj.WFSToolObj.reproduceAndRecordForAcousticPaths();
obj.WFSToolObj.calculateExperimentalAcousticPaths();
expAcPath = obj.WFSToolObj.expAcPathStruct;

%%
