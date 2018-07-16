%% Setup parameters

% Noise source variables
obj.amplitude = amplitude;
obj.amplitude(2) = -amplitude;
obj.phase = phase;

% Default values. They don't matter, but do not touch just in case.
obj.NSposition = [3.35 -0.2 0]; % Assumed real position
obj.frequency = 800;
obj.Fs = fs;

% Frequency filter
% obj.WFSToolObj.freqFilter = hTotal;
numFreqFilters = numel(freqFilters);
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

if ~exist('predefRoomDim', 'var')
    predefRoomDim = false;
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