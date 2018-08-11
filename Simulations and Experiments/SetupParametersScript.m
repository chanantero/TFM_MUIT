%% Setup parameters

% Noise source variables
if exist('amplitude', 'var')
    obj.amplitude = amplitude;
    obj.amplitude(2) = -amplitude;
else
    obj.amplitude = 1;
    obj.amplitude(2) = -1;
end

if exist('phase', 'var')
    obj.phase = phase;
else
    obj.phase = 0;
end

% Default values. They don't matter, but do not touch just in case.
obj.NSposition = [3.35 -0.2 0]; % Assumed real position
obj.frequency = 800;
obj.Fs = fs;

% Frequency filter
% obj.WFSToolObj.freqFilter = hTotal;
if exist('freqFilters', 'var')
    numFreqFilters = numel(freqFilters);
end
if ~exist('WFSfilterLength', 'var')
    predefWFSfilterLength = false;
else
    predefWFSfilterLength = true;
end
numFreqs = numel(freqs);

% Microphone positions
obj.microPos = recPositions;
numMicro = size(recPositions, 1);

% WFS options
if exist('frequencyCorrection', 'var')
    obj.WFSToolObj.frequencyCorrection = frequencyCorrection; % Very important! We want to see what happens without correction
else
    obj.WFSToolObj.frequencyCorrection = true;
end

if exist('attenuationType', 'var')
    obj.WFSToolObj.attenuationType = attenuationType;
else
    obj.WFSToolObj.attenuationType = 'Ruben';
end

% NS positions
numNSpos = size(NSpositions, 1);

% Room characteristics and impulse responses
Beta = beta(:) * [1 1 1 1 1 1];
numReverbTime = length(beta);

% Simulation options
if exist('automaticLengthModification', 'var')
    obj.WFSToolObj.automaticLengthModification = automaticLengthModification;
else
    obj.WFSToolObj.automaticLengthModification = true;
end

if ~exist('predefSignals', 'var')
    predefSignals = false;
end

if ~exist('roomDim', 'var')
    predefRoomDim = false;
else
    predefRoomDim = true;
end

if ~exist('predefNumSampIR', 'var')
    predefNumSampIR = false;
end
if ~exist('saveSignals', 'var')
    saveSignals = false;
end

if ~exist('progressBarActive', 'var')
    progressBarActive = true;
end