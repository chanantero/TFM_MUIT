%% Experiment 3
% Reproduction in GTAC

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

obj.setAcousticPaths('NS', 'theoretical');
%% Manual configuration of devices
% Select adequate audio driver and device for reproduction and recording

% Set number of microphones
numMicro = 1;
obj.WFSToolObj.setNumReceivers(numMicro);

% Set noise source channel through the GUI or use next line
NSchan = 1;
obj.WFSToolObj.noiseSourceChannelMapping(1) = NSchan;

%% Test that loudspeakers reproduce in an adequate way
numWFS = 96;
coefMat = zeros(numWFS, numWFS, 8);
semiToneShift = [-9 -7 -5 -4 -2 0 2 3];
freqLA = 440;
freqs = freqLA .* 2.^(semiToneShift/12);
prog = repmat([(1:7), (8:-1:2)], 1, ceil(numWFS/14));
indFreq = prog(1:numWFS);
indChan = 1:numWFS;
indPulse = 1:numWFS;
inds = sub2ind([numWFS, numWFS, 8], indPulse, indChan, indFreq);
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

%% Test for experimental frequency response
obj.WFSToolObj.reproduceAndRecordForAcousticPaths();
obj.WFSToolObj.calculateExperimentalAcousticPaths();
expAcPath = obj.WFSToolObj.expAcPathStruct;
obj.setLoudspeakerAcousticPath(expAcPath);

%% Calculate the coefficients for the cancellation
% Calculate coefficients for the cancellation with the experimental
% acoustic paths, but also with the official acoustic paths.

% Calculate WFS coefficients
obj.cancelResults = [];

    % WFS coefficients and global theoretical correction
obj.cancel({'NoFilter'}, true, {'Theoretical'}, {'AllTogether'}, false);

    % WFS coefficients and global correction with GTAC official impulse
% responses
load([dataPathName, 'acousticPathsGTAC_440.mat'])
obj.microPos = microphonePositions;
acPathWFSarrayStruct.acousticPaths = acousticPath;
acPathWFSarrayStruct.frequencies = 440;
obj.setAcousticPaths('NS', 'theoretical', 'WFS', acPathWFSarrayStruct);
obj.cancel({'NoFilter'}, true, {'Current'}, {'AllTogether'}, false);

    % WFS coefficients and global correction with experimental impulseResponses
obj.setLoudspeakerAcousticPath(expAcPath);
obj.cancel({'NoFilter'}, true, {'Current'}, {'AllTogether'}, false);
obj.setLoudspeakerAcousticPath(expAcPath);

% Simulate with experimental acoustic paths to see what we should get in
% the experimental case
obj.setLoudspeakerAcousticPath(expAcPath);
for k = 1:obj.numCancellationAttempts
    obj.WFScoef = obj.cancelResults(k).WFScoef;
    obj.WFSToolObj.prepareSimulation();
    obj.WFSToolObj.simulate();
    obj.cancelResults(k).recCoef = obj.microCoef;
    obj.cancelResults(k).recNScoef = obj.microCoefNS;
    obj.cancelResults(k).recWFScoef = obj.microCoefWFS;
end

%% Experimental checking
% Reproduce to check if the calculated predictions are veryfied by the
% experimental results
obj.experimentalChecking();
s = obj.cancelResults;

% Visualize 2D map
    % Format structure
for p = 1:numel(s)
    s(p).NScoef = [s(p).NSRcoef; s(p).NSVcoef];
    s(p).NSposition = [s(p).NSRposition; s(p).NSVposition];
end
    % Create simulationViewer object
objVis = simulationViewer(obj.ax, s);

