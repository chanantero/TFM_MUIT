% Create two noise sources. One of them will be virtual. It will be the
% assumed real source to be cancelled. The other will be real and its only
% purpose is to reproduce a signal.
% The real source will remain always the same. The virtual one will change
% in amplitude, phase and position.

%% Preamble
globalPath = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Matlab\';
paths = genpath(globalPath);
addpath(paths);

dataPathName = [globalPath, 'Data\'];
ID = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

%% System parameters

% User can change this:
numLoudspeakers = 96;

realPosition = [3.35 -0.2 0]; % Assumed real position
amplitude = 0.5;
phase = 0;
frequency = 440;

% Receivers
receiverPositions = [1.5 3 0; 1.5 1.5 0];
numReceivers = size(receiverPositions, 1);

% Predefined acoustic paths.
theoricWFSAcPath = true;
theoricNoiseSourceAcPath = true;
% acPathWFSarray;
% acPathNoiseSources;

%% Set Up System
if exist('obj', 'var') == 0 || ~isvalid(obj) || ~isvalid(obj.ax)
    obj = WFSToolSimple;
end

% Set noise source variables
obj.setNumNoiseSources(2);
obj.noiseSourceChannelMapping = [1; 0];
obj.amplitude = [amplitude; amplitude];
obj.phase = [phase; phase];
obj.frequency = [frequency; frequency];
obj.noiseSourcePosition = [realPosition; realPosition];
obj.setVirtual([false; true]);
obj.setReal([true; false]);

% Set WFS array variables
obj.setNumWFSarraySources(numLoudspeakers);

% Set receiver variables
obj.setNumReceivers(numReceivers);
obj.receiverPosition = receiverPositions;

% Set acoustic path variables
if theoricWFSAcPath
    obj.theoricWFSacousticPath();
else
    obj.WFSarrayAcPathStruct.acousticPaths = acPathWFSarray;
end

if theoricNoiseSourceAcPath
    obj.theoricNoiseSourceAcousticPath();
else
    obj.noiseSourceAcPathStruct.acousticPaths = acPathNoiseSources;
end

obj.updateReprodPanelBasedOnVariables();
obj.updateRecordPanelBasedOnVariables();

%% Experimental acoustic path calculation

% Calculate experimental acoustic path
obj.reproduceAndRecordForAcousticPaths();
obj.calculateExperimentalAcousticPaths();

% Retrieve and save information
sBase = obj.exportInformation();
FileName = ['Session_', ID, 'calibration', '.mat'];
save([dataPathName, FileName], 'sBase');

%% Set experimental acoustic path
ID = '2017-09-19_18-18-47';
fileName = ['Session_', ID, 'calibration', '.mat'];
s = load([dataPathName, fileName], '-mat' , 'sBase');
sBase = s.sBase;
expAcPath = sBase.Experiment.acPathStruct;
obj.setLoudspeakerAcousticPath(expAcPath);
obj.noiseSourceAcPathStruct.acousticPaths(:, 2) = obj.noiseSourceAcPathStruct.acousticPaths(:, 1); % The acoustic path of the virtual noise source is equal to the real one, for optimization purposes

%% WFS cancellation

% Optimization options
sourceFilter = {'Loudspeakers', 'Loudspeakers'}; % It makes no sense to optimize a less real scenario, so use loudspeakers filter
maxAbsValCons = [true, true]; % We always want to be realistic about real constraints
acousticPathType = {'Current', 'Current'}; % It makes no sense to optimize with a theoric acoustic path because it depends on the parameters of the noise source, and those parameters are actually unknown. Besides, the acousic path of the loudspeakers is only known in the places where microphones have been placed.
grouping = {'Independent', 'Independent'}; 
zerosFixed = [true, false];
N = numel(sourceFilter);

loudsCoeff = zeros(obj.numLoudspeakers, obj.numNoiseSources, N);
simulField = zeros(obj.numReceivers, obj.numNoiseSources, N);
for k = 1:N
% WFS cancellation
obj.WFScalculation('SourceFilter', sourceFilter{k}, 'AcousticPath', acousticPathType{k}, 'Grouping', grouping{k}, 'maxAbsoluteValueConstraint', maxAbsValCons(k), 'zerosFixed', zerosFixed(k));

% WFS array coefficients
loudsCoeff(:, :, k) = obj.loudspeakerCoefficient;

% Simulate
obj.simulate();

% Cancellation level (noise source 1 is real, noise source 2 is virtual)
simulField(:, :, k) = obj.simulField;
end

noiseSourceSimulField = permute(simulField(:, 1, :), [1 3 2]); % (numReceivers x N)
totalSimulField = permute(sum(simulField, 2), [1 3 2]); % (numReceivers x N)
simulRelField = totalSimulField./noiseSourceSimulField;
simulCancel = 20*log10(abs(simulRelField));

%% Experimental checking of correspondence

% A) Reproduce
% A.1) Only noise
obj.setVirtual([false; false]); % No virtual sources
obj.setReal([true; false]); % Only real source
obj.updateReprodPanelBasedOnVariables();

obj.WFScalculation(); % Update coefficents of loudspeakers. Only the channel of the real loudspeaker should have coefficient different than 0
obj.reproduceAndRecord('main', 'soundTime', 2); % Simple reproduction of one pulse of 2 seconds

% Retrieve and save information
sExpOnlyNoise = obj.getExperimentalResultVariables();
yPulseCoefMat_OnlyNoise = signal2pulseCoefficientMatrix(sExpOnlyNoise.pulseLimits, frequency, ones(1, 1), sExpOnlyNoise.recordedSignal, sExpOnlyNoise.sampleRate);

FileName = ['Session_', ID, 'onlyNoise', '.mat'];
save([dataPathName, FileName], 'sExpOnlyNoise', 'yPulseCoefMat_OnlyNoise');

% Go back to previous configuration
obj.setVirtual([false; true]);
obj.setReal([true; false]);
obj.updateReprodPanelBasedOnVariables();
obj.WFScalculation();

% A.2) Noise + WFS array
% Create signal
pulseDuration = 1;
silenceDuration = 1;
startPulse = (0:N - 1)'*(pulseDuration + silenceDuration);
endPulse = startPulse + pulseDuration;
pulseLim = [startPulse, endPulse];
pulseCoefMat = permute(loudsCoeff, [3 1 2]);
SampleRate = 44100;
signalFunc = @(startSample, endSample) pulseCoefMat2signal(pulseCoefMat, pulseLim, obj.frequency, SampleRate, startSample, endSample, 'type_pulseLimits', 'time');

% Reproduce
obj.reproduceSignalFunction(signalFunc, SampleRate, obj.loudspeakerChannelMapping);

obj.pulseCoeffMat = pulseCoefMat;
obj.pulseLimits = pulseLim;
obj.reprodFrequencies = obj.frequency;

% Retrieve and save information
sExp = obj.getExperimentalResultVariables();
yPulseCoefMat = signal2pulseCoefficientMatrix(sExp.pulseLimits, sExp.frequencies(1), sum(sExp.pulseCoefMat, 3), sExp.recordedSignal, sExp.sampleRate);

FileName = ['Session_', ID, 'cancellation', '.mat'];
save([dataPathName, FileName], 'sExp', 'yPulseCoefMat');

% B) Check how similar the results are with the simulation
% yPulseCoefMat % (N x numReceivers)
% simulField % (numReceivers x numNoiseSources x N)
noiseSourceExpField = repmat(yPulseCoefMat_OnlyNoise.', [1, N]);
totalExpField = permute(yPulseCoefMat, [2 1]); % (numReivers x N)

ratio = totalExpField./totalSimulField;
maxAbsRat = max(abs(ratio));
maxAngleRat = max(angle(ratio)) - min(angle(ratio));

expRelField = totalExpField./noiseSourceExpField;
expCancel = 20*log10(abs(expRelField));


%% Optimization of noise source theoric parameters
% Find the theoric parameters for the virtual noise source that, applying
% WFS cancellation with theoric acoustic path and unified optimization of loudspeakers, minimize the
% magnitude of the resulting field using the experimental acoustic path

% Optimize
% The objective function accepts parameters as input arguments. It
% returns the absolute value of the resulting field
objectiveFunction = @(parameters) sum(abs(sum(noiseSourceParam2Field(obj, parameters), 2)).^2);
x0 = [realPosition, amplitude, phase]; % Initial value of parameters
[xOpt, fVal] = fminunc(objectiveFunction, x0); % Optimize

% Set the parameters
optPosition = xOpt(1:3);
optAmplitude = xOpt(4);
optPhase = xOpt(5);
% optPosition = x0(1:3);
% optAmplitude = x0(4);
% optPhase = x0(5);
obj.amplitude(2) = optAmplitude;
obj.phase(2) = optPhase;
obj.noiseSourcePosition(2, :) = optPosition;
obj.updateReprodPanelBasedOnVariables();

% Set WFS array coefficients according to the optimized parameters. Respect
% the generated theoric coefficients, don't ajust them independently or
% something, just scale them to fit the theoric acoustic paths.
obj.WFScalculation(...
    'SourceFilter', 'Loudspeakers',...
    'AcousticPath', 'Theoric',...
    'Grouping', 'AllTogether',...
    'maxAbsoluteValueConstraint', true);

% Simulate with the experimental acoustic path to see what is the resulting
% field with these optimized theoric parameters of the virtual source
obj.simulate();

% Calculate cancellation
noiseSourceField = obj.simulField(:, 1); % (numReceivers x 1)
totalField = sum(obj.simulField, 2); % (numReceivers x 1)
relField = totalField./noiseSourceField;
cancel = 20*log10(abs(relField));

%% Experimental checking of correspondence
% Reproduction of the opmized virtual noise source theoric parameters
% Variables that should be ready at the beginning of this section
% pulseCoefMat. (N x obj.numLoudspeakers x obj.numNoiseSources)
pulseCoefMat = permute(obj.loudspeakerCoefficient, [3 1 2]);


% A) Reproduction
% A.2) Noise + WFS array
% Create signal
pulseDuration = 1;
silenceDuration = 1;
numPulses = size(pulseCoefMat, 1);
startPulse = (0:numPulses - 1)'*(pulseDuration + silenceDuration);
endPulse = startPulse + pulseDuration;
pulseLim = [startPulse, endPulse];
SampleRate = 44100;
signalFunc = @(startSample, endSample) pulseCoefMat2signal(pulseCoefMat, pulseLim, obj.frequency, SampleRate, startSample, endSample, 'type_pulseLimits', 'time');

% Reproduce
obj.reproduceSignalFunction(signalFunc, SampleRate, obj.loudspeakerChannelMapping);

obj.pulseCoeffMat = pulseCoefMat;
obj.pulseLimits = pulseLim;
obj.reprodFrequencies = obj.frequency;

% Retrieve and save information
sExp = obj.getExperimentalResultVariables();
yPulseCoefMat = signal2pulseCoefficientMatrix(sExp.pulseLimits, sExp.frequencies(1), sum(sExp.pulseCoefMat, 3), sExp.recordedSignal, sExp.sampleRate);

FileName = ['Session_', ID, 'optimizedParam', '.mat'];
save([dataPathName, FileName], 'sExp', 'yPulseCoefMat');

% B) Check how similar the results are with the simulation
% yPulseCoefMat % (numPulses x numReceivers)
% simulField % (numReceivers x numNoiseSources x numPulses)
noiseSourceExpField = repmat(yPulseCoefMat_OnlyNoise.', [1, numPulses]);
totalExpField = permute(yPulseCoefMat, [2 1]); % (numReivers x numPulses)

ratio = totalExpField./totalSimulField;
maxAbsRat = max(abs(ratio));
maxAngleRat = max(angle(ratio)) - min(angle(ratio));

expRelField = totalExpField./noiseSourceExpField;
expCancel = 20*log10(abs(expRelField));


