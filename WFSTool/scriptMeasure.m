% Create two noise sources. One of them will be virtual. It will be the
% assumed real source to be cancelled. The other will be real and its only
% purpose is to reproduce a signal.
% The real source will remain always the same. The virtual one will change
% in amplitude, phase and position.
obj = WFSToolSimple;

numLoudspeakers = 96;

realPosition = [0.6 0 0]; % Assumed real position

minXPos = realPosition(1); maxXPos = realPosition(1); numXPoints = 1;
minYPos = realPosition(2); maxYPos = realPosition(2); numYPoints = 1;
minZPos = 0; maxZPos = 0; numZPoints = 1;
minAmplitude = 0; maxAmplitude = 1; numAmplitudePoints = 5;
numPhasePoints = 1;

amplitude = 0.5;
phase = 0;
frequency = 600;

% Set scenario
obj.changeScenario(numLoudspeakers);

obj.setNumNoiseSources(2);
obj.noiseSourceChannelMapping = [1; 0];

obj.amplitude = [amplitude; amplitude];
obj.phase = [phase; phase];
obj.frequency = [frequency; frequency];
obj.noiseSourcePosition = [realPosition; realPosition];

obj.updateReprodPanelBasedOnVariables();

% For the current configuration. Get the experimental acoustic path
obj.reproduceAndRecordForAcousticPaths();
obj.calculateExperimentalAcousticPaths();
sBase = obj.exportInformation();

PathName = 'Data\';
FileName = sprintf('CalibrationLab_31-08-2017.mat');
save([PathName, FileName], 'sBase');

% Test only the real noise source.
obj.setVirtual([false; false]);
obj.setReal([true; false]);
obj.updateReprodPanelBasedOnVariables();
obj.WFScalculation();
obj.reproduceAndRecord('main');
sExpOnlyNoise = obj.getExperimentalResultVariables();
yPulseCoefMat_OnlyNoise = signal2pulseCoefficientMatrix(frequency, ones(1, 1), sExpOnlyNoise.pulseLimits, sExpOnlyNoise.recordedSignal, sExpOnlyNoise.sampleRate);

obj.setVirtual([false; true]);
obj.setReal([true; false]);
obj.updateReprodPanelBasedOnVariables();
obj.WFScalculation();

% Test the multiple cases
xVec = linspace(minXPos, maxXPos, numXPoints);
yVec = linspace(minYPos, maxYPos, numYPoints);
zVec = linspace(minZPos, maxZPos, numZPoints);
amplitudeVec = linspace(0, 1, numAmplitudePoints);
phaseVec = 2*pi/numPhasePoints*(0:numPhasePoints - 1);
[X, Y, Z, A, P] = ndgrid(xVec, yVec, zVec, amplitudeVec, phaseVec);
sizeMat = size(X);
numPointsMat = numel(X);

pulseCoefMat = zeros(numPointsMat, numLoudspeakers, 2);
for p = 1:numPointsMat
    % Set virtual position
    virtPos = [X(p), Y(p), Z(p)];
    obj.noiseSourcePosition = [realPosition; virtPos];
    
    % Set amplitude and phase
    obj.amplitude(2) = A(p);
    obj.phase(2) = P(p);
    
    % Apply WFS calculation
    obj.WFScalculation();
    
    % Get the coefficients and save them in the pulse coefficient matrix
    pulseCoefMat(p, :, :) = permute(obj.loudspeakerCoefficient, [3, 1, 2]);
end
    
% Create signal
pulseDuration = 1;
silenceDuration = 1;
startPulse = (0:numPointsMat - 1)'*(pulseDuration + silenceDuration);
endPulse = startPulse + pulseDuration;
pulseLim = [startPulse, endPulse];
SampleRate = 44100;
signalFunc = @(startSample, endSample) pulseCoefMat2signal(obj.frequency, pulseCoefMat, floor(pulseLim*SampleRate), SampleRate, startSample, endSample, 'sample');

% Reproduce
obj.reproduceSignalFunction(signalFunc, SampleRate, obj.loudspeakerChannelMapping);

obj.pulseCoeffMat = pulseCoefMat;
obj.pulseLimits = pulseLim;
obj.reprodFrequencies = obj.frequency;

% Retrieve and save information
PathName = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Matlab\Data\';
FileName = sprintf('Example.mat');

sExp = obj.getExperimentalResultVariables();
save([PathName, FileName], 'sExp');

%% Analyse results
yPulseCoefMat = signal2pulseCoefficientMatrix(sExp.pulseLimits, sExp.frequencies(1), sum(sExp.pulseCoefMat, 3), sExp.recordedSignal, sExp.sampleRate);

% Simulate pulse coefficients based on experimental acoustic path
expAcPath = tuneAcousticPaths(sExp.acPathStruct.acousticPaths, sExp.acPathStruct.frequencies, frequency);
numReceivers = obj.numReceivers;
yPulseCoefMat_pseudoExp = zeros(numPointsMat, numReceivers, 2);
for p = 1:numPointsMat
    yPulseCoefMat_pseudoExp(p, :, :) = sum(expAcPath .* repmat(permute(pulseCoefMat(p, :, :), [2, 1, 3]), [1, numReceivers, 1]), 1);    
end
yPulseCoefMat_pseudoExp = sum(yPulseCoefMat_pseudoExp, 3);

% Compare real pulse coefficients with the pseudo-experimental ones

% Use predefined acoustic path
acPathWFSarray;
acPathNoiseSources;
receiverPositions;

acPath = acPathWFSarray;
[~, indReal] = ismember(obj.noiseSourceChannelMapping(obj.real), obj.loudspeakerChannelMapping);
acPath(:, indReal, :) = acPathNoiseSources;

yPulseCoefMat = zeros(numPointsMat, numReceivers, 2);
for f = 1:numFrequencies
    yPulseCoefMat(:, :, f) = pulseCoefMat(:, :, f)*acPath.';
end
yPulseCoefMat = sum(yPulseCoefMat, 3);

obj.simulObj.measurePoints = receiverPositions;
for p = 1:numPulses
    obj.simulObj.field = permute(yPulseCoefMat(p, :, :), [2 3 1]);
    obj.simulObj.draw('scatter');
end

