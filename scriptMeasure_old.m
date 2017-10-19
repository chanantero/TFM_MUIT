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

%% Parameters.

% User can change this:
numLoudspeakers = 96;

realPosition = [3.35 -0.2 0]; % Assumed real position
amplitude = 0.5;
phase = 0;
frequency = 440;

minXPos = 3.35; maxXPos = 3.35; numXPoints = 1;
minYPos = -0.2; maxYPos = -0.2; numYPoints = 1;
minZPos = 0; maxZPos = 0; numZPoints = 1;
minAmplitude = 0; maxAmplitude = 1; numAmplitudePoints = 11;
numPhasePoints = 1;

% Receivers
receiverPositions = [1.5 3 0; 1.5 1.5 0];
numReceivers = size(receiverPositions, 1);

% Predefined acoustic paths.
theoricWFSAcPath = true;
theoricNoiseSourceAcPath = true;
% acPathWFSarray;
% acPathNoiseSources;

%% Set scenario
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

%% Calculate useful variables for the different cases
% Calculate the multiple cases
xVec = linspace(minXPos, maxXPos, numXPoints);
yVec = linspace(minYPos, maxYPos, numYPoints);
zVec = linspace(minZPos, maxZPos, numZPoints);
amplitudeVec = linspace(minAmplitude, maxAmplitude, numAmplitudePoints);
phaseVec = 2*pi/numPhasePoints*(0:numPhasePoints - 1);

[X, Y, Z, A, P] = ndgrid(xVec, yVec, zVec, amplitudeVec, phaseVec);
numPointsMat = numel(X);

% For the visualization we want to transform the results into an array with
% each dimension separated. Specifically, we want to rearrange data into
% a multidimensional array:
% 1st dimension: x position of the virtual noise source
% 2nd dimension: y position of the virtual noise source
% 3rd dimension: z position of the virtual noise source
% 4th dimension: amplitude vector
% 5th dimension: phase vector
% 6th dimension: channels

% Use the function mergeAndPermute in inverse mode (third argument to true)
% The operation that we would have to do to transform that array into what
% we have is described by modVec:
modVec = {[1, 2, 3, 4, 5], 6};
% Size of the array
sizeObjective = [numXPoints, numYPoints, numZPoints, numAmplitudePoints, numPhasePoints, numReceivers];
% Label of each dimension
labels = {'x', 'y', 'z', 'Amplitude', 'Phase', 'Microphone'};

%% Simulation of cases
% Calculate the received field only with the real source, without
% cancellation
obj.setVirtual([false; false]); % No virtual sources
obj.WFScalculation();
obj.simulate();
onlySourceSimulField = obj.simulField(:, 1);
obj.setVirtual([false; true]); % No virtual sources


pulseCoefMat = zeros(numPointsMat, obj.numLoudspeakers, 2);
simulatedField = zeros(obj.numReceivers, obj.numNoiseSources, numPointsMat);
for p = 1:numPointsMat
    fprintf('Point %d\n', p);
    
    % Set properties
        % Set virtual position
        virtPos = [X(p), Y(p), Z(p)];
        obj.noiseSourcePosition = [realPosition; virtPos];
        obj.theoricNoiseSourceAcousticPath();

        % Set amplitude and phase
        obj.amplitude(2) = A(p);
        obj.phase(2) = P(p);

        % Apply WFS calculation
        obj.WFScalculation();

    % Simulate the received signals
    obj.simulate();
    simulatedField(:, :, p) = obj.simulField;

    % Get the coefficients and save them in the pulse coefficient matrix
    pulseCoefMat(p, :, :) = permute(obj.loudspeakerCoefficient, [3, 1, 2]);
    
end
simulatedField = sum(permute(simulatedField, [3, 1, 2]), 3);
simulatedFieldRelative = simulatedField./repmat(onlySourceSimulField.', [numPointsMat, 1]);
simulatedFieldRelative_DB = 20*log10(abs(simulatedFieldRelative));

%% Experiment

% 1) Experimental acoustic path
obj.reproduceAndRecordForAcousticPaths();
obj.calculateExperimentalAcousticPaths();

% Retrieve and save information
sBase = obj.exportInformation();
FileName = ['Session_', ID, 'calibration', '.mat'];
save([dataPathName, FileName], 'sBase');

% 2) Only real noise source
obj.setVirtual([false; false]); % No virtual sources
obj.setReal([true; false]); % Only real source
obj.updateReprodPanelBasedOnVariables();

obj.WFScalculation(); % Update coefficents of loudspeakers. Only the channel of the real loudspeaker should have coefficient different than 0
obj.reproduceAndRecord('main', 'soundTime', 2); % Simple reproduction of one pulse of 1 second

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

% 3) Reproduce signal for the multiple cases

% Create signal
pulseDuration = 1;
silenceDuration = 1;
startPulse = (0:numPointsMat - 1)'*(pulseDuration + silenceDuration);
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

FileName = ['Session_', ID, 'cases', '.mat'];
save([dataPathName, FileName], 'sExp', 'yPulseCoefMat');

%% Visualize
visualizeScript;

%% Analyse Experimental Results
% Get experimental data
    % General information structure after calibration
    s = load('C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Matlab\Data\Session_2017-09-19_18-18-47calibration.mat', '-mat' , 'sBase');
    sBase = s.sBase;
    % Acoustic path
    expAcPath = sBase.Experiment.acPathStruct;
    % Only the noise source
    s = load('C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Matlab\Data\Session_2017-09-19_18-18-47onlyNoise.mat', '-mat' , 'sExpOnlyNoise', 'yPulseCoefMat_OnlyNoise');
    sExpOnlyNoise = s.sExpOnlyNoise;
    yPulseCoefMat_OnlyNoise = s.yPulseCoefMat_OnlyNoise; % (numPulses x numReceivers)
    % Cases
    s = load('C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Matlab\Data\Session_2017-09-19_18-18-47cases.mat', '-mat' , 'sExp', 'yPulseCoefMat');
    sExpCases = s.sExp;
    yPulseCoefMat_Cases= s.yPulseCoefMat;

    clear('s')

% Check that the experimental results are the same (approximatelly) than
% the ones simulated with the experimental acoustic path.

% Set the experimental acoustic paths
obj.setLoudspeakerAcousticPath(expAcPath);
obj.WFSarrayCoefficient(:) = 0; % Just in case
obj.noiseSourceCoefficient_complete(:) = 0; % Just in case

% Only the noise source simulation
pulseCoefMat = sExpOnlyNoise.pulseCoefMat;
loudsChanMap = sExpOnlyNoise.channelMapping;

numPointsMat = size(pulseCoefMat, 1);
simulatedField = zeros(obj.numReceivers, obj.numNoiseSources, numPointsMat);
for p = 1:numPointsMat
    loudsCoef = permute(pulseCoefMat(p, :, :), [2, 3, 1]);
    obj.setLoudspeakerCoefficient(loudsCoef, loudsChanMap);
    obj.simulate();
    simulatedField(:, :, p) = obj.simulField;
end

yPulseCoefMat_OnlyNoise_Simulated = sum(permute(simulatedField, [3, 1, 2]), 3);

% Cases simulation
pulseCoefMat = sExpCases.pulseCoefMat;
loudsChanMap = sExpCases.channelMapping;

numPointsMat = size(pulseCoefMat, 1);
simulatedField = zeros(obj.numReceivers, obj.numNoiseSources, numPointsMat);
for p = 1:numPointsMat
    loudsCoef = permute(pulseCoefMat(p, :, :), [2, 3, 1]);
    obj.setLoudspeakerCoefficient(loudsCoef, loudsChanMap);
    obj.simulate();
    simulatedField(:, :, p) = obj.simulField;
end

yPulseCoefMat_Cases_Simulated = sum(permute(simulatedField, [3, 1, 2]), 3);

% Compare the simulated fields with the real ones
rat_OnlyNoise = yPulseCoefMat_OnlyNoise_Simulated./yPulseCoefMat_OnlyNoise;
rat_Cases = yPulseCoefMat_Cases_Simulated./yPulseCoefMat_Cases;

% It seems it's approximatelly the same. HURRA!!

% Visualize cancellation in dB of the experimental field
relativeField = yPulseCoefMat_Cases./repmat(yPulseCoefMat_OnlyNoise(1, :), [numPointsMat, 1]);
relativeField_dB = 20*log10(abs(relativeField));
data = relativeField_dB;
visualizeScript;

%% Least squares for cancellation
% Cancel the field using the experimental acoustic path and the least
% squares method
    % Set the virtual noise coefficient to the same value as the real one, as
    % the cancellation happens with each component/source/frequency
    % independently
    obj.noiseSourceCoefficient(2) = obj.noiseSourceCoefficient(1);
    
    % Make sure the current acoustic path is the experimental one
    obj.setLoudspeakerAcousticPath(expAcPath);

    % Apply least squares method
    obj.WFScalculation('SourceFilter', 'Loudspeakers', 'AcousticPath', 'Current', 'Grouping', 'Independent');

% Simulate to see if, indeed, the field has been cancelled
obj.simulate();
sum(obj.simulField, 2) % Sum along noise source components

%% Least squares for cancellation with the official impulse responses of the GTAC
% acousticPath = importImpulseResponseGTAC(frequency);
load([dataPathName, 'acousticPathsGTAC_440.mat'])

% Set receiver properties
obj.setNumReceivers(size(microphonePositions, 1)); % Doesn't matter the position, although we know them
obj.receiverPosition = microphonePositions;
obj.updateRecordPanelBasedOnVariables();

% Acoustic Path Assignation
obj.theoricNoiseSourceAcousticPath();
obj.WFSarrayAcPathStruct.acousticPaths = acousticPath;

% Optimization options
sourceFilter = {'Loudspeakers', 'Loudspeakers', 'Loudspeakers'};
acousticPathType = {'Current', 'Current', 'Current'};
grouping = {'Independent', 'AllTogether', 'AllTogether'};
maxAbsValCons = [false, false, true];
zerosFixed = [true, true, true];
N = numel(sourceFilter);

simulField = zeros(obj.numReceivers, obj.numNoiseSources, N);
for k = 1:N
% WFS cancellation
obj.WFScalculation('SourceFilter', sourceFilter{k}, 'AcousticPath', acousticPathType{k}, 'Grouping', grouping{k}, 'maxAbsoluteValueConstraint', maxAbsValCons(k), 'zerosFixed', zerosFixed(k));

% Simulate
obj.simulate();

% Cancellation level (noise source 1 is real, noise source 2 is virtual)
simulField(:, :, k) = obj.simulField;
end

noiseSourceField = permute(simulField(:, 1, :), [1 3 2]); % (numReceivers x N)
totalField = permute(sum(simulField, 2), [1 3 2]); % (numReceivers x N)
relField = totalField./noiseSourceField;
cancel_dB = 20*log10(abs(relField));

ax = axes(figure);
plot(ax, cancel_dB, 'Marker', '.')
plot(ax, [abs(totalField), abs(noiseSourceField)])



