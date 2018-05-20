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
obj.WFSToolObj.filterWFS_length = WFSfilterLength;

% Microphone positions
obj.microPos = recPositions;
numMicro = size(recPositions, 1);

% NS positions
numNSpos = size(NSpositions, 1);

% Room characteristics and impulse responses
Beta = beta(:) * [1 1 1 1 1 1];
numReverbTime = length(beta);

% WFS options
obj.WFSToolObj.frequencyCorrection = frequencyCorrection; % Very important! We want to see what happens without correction

% Simulation options
obj.WFSToolObj.automaticLengthModification = automaticLengthModification;