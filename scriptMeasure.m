% Create two noise sources. One of them will be virtual. It will be the
% assumed real source to be cancelled. The other will be real and its only
% purpose is to reproduce a signal.
% The real source will remain always the same. The virtual one will change
% in amplitude, phase and position.

%% Parameters.
paths = genpath('C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Matlab\');
addpath(paths);

% User can change this:
numLoudspeakers = 2;

realPosition = [1 -0.2 0]; % Assumed real position
amplitude = 0.5;
phase = 0;
frequency = 600;

minXPos = 0.8; maxXPos = 1.2; numXPoints = 10;
minYPos = -0.2; maxYPos = -0.2; numYPoints = 1;
minZPos = 0; maxZPos = 0; numZPoints = 1;
minAmplitude = 0.5; maxAmplitude = 0.5; numAmplitudePoints = 1;
numPhasePoints = 1;

% Receivers
receiverPositions = [1.5 3 0; 1.5 1.5 0];

% Predefined acoustic paths.
theoricWFSAcPath = true;
theoricNoiseSourceAcPath = true;
% acPathWFSarray;
% acPathNoiseSources;

PathName = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Matlab\Data\';
ID = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

%% Processing

%% A. Set scenario
if exist('obj', 'var') == 0 || ~isvalid(obj) || ~isvalid(obj.ax)
    obj = WFSToolSimple;
end

obj.setNumWFSarraySources(numLoudspeakers);
obj.setNumReceivers(size(receiverPositions, 1));

obj.setNumNoiseSources(2);
obj.noiseSourceChannelMapping = [1; 0];
obj.amplitude = [amplitude; amplitude];
obj.phase = [phase; phase];
obj.frequency = [frequency; frequency];
obj.noiseSourcePosition = [realPosition; realPosition];
obj.setVirtual([false; true]);
obj.setReal([true; false]);
obj.setNumReceivers(size(receiverPositions, 1));
obj.receiverPosition = receiverPositions;

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
%% B. Preparte the variables for the multiple cases
% Calculate the received field only with the real source, without
% cancellation
obj.setVirtual([false; false]); % No virtual sources
obj.WFScalculation();
obj.simulate();
onlySourceSimulField = obj.simulField(:, 1);
obj.setVirtual([false; true]); % No virtual sources

% Calculate the multiple cases
xVec = linspace(minXPos, maxXPos, numXPoints);
yVec = linspace(minYPos, maxYPos, numYPoints);
zVec = linspace(minZPos, maxZPos, numZPoints);
amplitudeVec = linspace(minAmplitude, maxAmplitude, numAmplitudePoints);
phaseVec = 2*pi/numPhasePoints*(0:numPhasePoints - 1);
[X, Y, Z, A, P] = ndgrid(xVec, yVec, zVec, amplitudeVec, phaseVec);
sizeMat = size(X);
numPointsMat = numel(X);

pulseCoefMat = zeros(numPointsMat, obj.numLoudspeakers, 2);
simulatedField = zeros(obj.numReceivers, obj.numNoiseSources, numPointsMat);
for p = 1:numPointsMat
    fprintf('Point %d\n', p);
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

%% C. Experiment

% C.1 Experimental acoustic path
obj.reproduceAndRecordForAcousticPaths();
obj.calculateExperimentalAcousticPaths();

% Retrieve and save information
sBase = obj.exportInformation();
FileName = ['Session_', ID, 'calibration', '.mat'];
save([PathName, FileName], 'sBase');

% C.2 Only real noise source
obj.setVirtual([false; false]); % No virtual sources
obj.setReal([true; false]); % Only real source
obj.updateReprodPanelBasedOnVariables();

obj.WFScalculation(); % Update coefficents of loudspeakers. Only the channel of the real loudspeaker should have coefficient different than 0
obj.reproduceAndRecord('main', 'soundTime', 1); % Simple reproduction of one pulse of 1 second

% Retrieve and save information
sExpOnlyNoise = obj.getExperimentalResultVariables();
yPulseCoefMat_OnlyNoise = signal2pulseCoefficientMatrix(sExpOnlyNoise.pulseLimits, frequency, ones(1, 1), sExpOnlyNoise.recordedSignal, sExpOnlyNoise.sampleRate);

FileName = ['Session_', ID, 'onlyNoise', '.mat'];
save([PathName, FileName], 'sExpOnlyNoise', 'yPulseCoefMat_OnlyNoise');

% Go back to previous configuration
obj.setVirtual([false; true]);
obj.setReal([true; false]);
obj.updateReprodPanelBasedOnVariables();
obj.WFScalculation();

% C.3 Reproduce signal for the multiple cases

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
save([PathName, FileName], 'sExp', 'yPulseCoefMat');


%%
% Simulate pulse coefficients based on experimental acoustic path


% Compare real pulse coefficients with the pseudo-experimental ones


%% Visualize

% Rearrange data into multidimensional array
% 1st dimension: x position of the virtual noise source
% 2nd dimension: y position of the virtual noise source
% 3rd dimension: z position of the virtual noise source
% 4th dimension: amplitude vector
% 5th dimension: phase vector
% 6th dimension: channels

% The operation to obtain simulatedField from the described array is:
modVec = {[1, 2, 3, 4, 5], 6};
sizeFormatted = [numXPoints, numYPoints, numZPoints, numAmplitudePoints, numPhasePoints, obj.numReceivers];
simulFieldFormatted = mergeAndPermute(simulatedFieldRelative, modVec, true, sizeFormatted);

visualObj = animation({xVec, yVec, zVec, amplitudeVec, phaseVec, 1:obj.numReceivers},...
    {20*log10(abs(simulFieldFormatted))}, {'x', 'y', 'z', 'Amplitude', 'Phase', 'Microphone'}, {'Cancellation (dB)'}, [], []);


%% Analyse
sExp = sBase.Experiment;
[expAcPath, freq, yPulseCoefMat] = getAcousticPath( sExp.pulseLimits, sExp.frequencies, sExp.pulseCoefMat, sExp.recordedSignal, sExp.sampleRate);

bar(abs(expAcPath)')
bar(rad2deg(angle(expAcPath))')
