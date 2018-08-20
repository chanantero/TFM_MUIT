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


freqs = 0;

% As the audio driver only allows 96 channels and we need the reproduction
% of the noise source signal to be synchronized, it is necessary to
% dedicate one of the channels to the noise source. Of course, the corresponding
% array loudspeaker will be in silence.
% Set noise source channel through the GUI or use next line
NSchan = 1;
obj.WFSToolObj.noiseSourceChannelMapping(1) = NSchan;

obj.domain = 'time';

zPos = 1.65;
    WFSarrayOffset = [0.46 2.21 zPos]; % [x, y, z] coordinates. Useful for generating acoustic path IR.
    roomDim = [4.48, 9.13, 2.64];
    fs = 44100/4;
    c = 340;
    
    extRectXmin = min(obj.WFSposition(:, 1));
    extRectXmax = max(obj.WFSposition(:, 1));
    extRectYmin = min(obj.WFSposition(:, 2));
    extRectYmax = max(obj.WFSposition(:, 2));
    centerX = (extRectXmax + extRectXmin)/2;
    centerY = (extRectYmax + extRectYmin)/2;
    recPositions = [centerX, centerY, 0; centerX, centerY + 1, 0];
    
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
    NSpositions = obj.NSRposition;
    
     % WFS options
    frequencyCorrection = true;
    attenuationType = 'Ruben';
    
    % Simulation options
    timeDomainActive = true;
    fakeTimeProcessing = false;
    frequencyDomainActive = true;
    automaticLengthModification = false;

    
    % Define chirp signal
durSign = 4; % Duration of tone for time processing
sampleRate = fs;
t = (0:ceil(durSign*sampleRate)-1)/sampleRate;
NSsignal = chirp(t, 20, durSign, 950);
prefixDuration = 2;
prefixNumSamples = floor(prefixDuration*sampleRate) + 1;
NSsignal = [zeros(1, prefixNumSamples), NSsignal, zeros(1, prefixNumSamples )];
t = (0:length(NSsignal)-1)/sampleRate;
predefSignals = true;
saveSignals = true;

% Frequency filters
magnFiltOrder = 2^10;
hilbertFiltOrder = 2^12;
[freqFilter, freqFiltDelay] = getFrequencyFilter(magnFiltOrder, hilbertFiltOrder, fs);
freqFilters = {freqFilter};
freqFiltDelays = freqFiltDelay;
    
SetupParametersScript
AcousticPathCalculationScript
simulationScript
    
ax = axes(figure);
t = (0:size(recNS_signals, 2)-1)/sampleRate;
plot(ax, t, recNS_signals)

%%

obj = SimulationController;

% Constants
c = 340; % Sound velocity (m/s)
fs = 44100; % Sample frequency (samples/s)
WFSarrayOffset = [0.5, 1.5, 1.5]; % [x, y, z] coordinates. Useful for generating acoustic path IR.

% Noise source coefficient
amplitude = 1;
phase = 0;

% Filter variables for the time WFS filter.
WFSfilterLength = 22050;
% Creation of frequency filters with different orders.
magnFiltOrder = 2.^(12);
hilbertFiltOrder = 2.^(12);
numFreqFilters = length(magnFiltOrder);

freqFilters = cell(numFreqFilters, 1);
freqFiltDelays = zeros(numFreqFilters, 1);
for k = 1:numFreqFilters
    [freqFilter, delay] = getFrequencyFilter( magnFiltOrder(k), hilbertFiltOrder(k), fs );    
    freqFilters{k} = freqFilter;
    freqFiltDelays(k) = delay;
end

% Microphone positions
% Rectangular grid
marginRatio = 0.3;
numPointsX = 2;
numPoinstY = 2;
extRectXmin = min(obj.WFSposition(:, 1));
extRectXmax = max(obj.WFSposition(:, 1));
extRectYmin = min(obj.WFSposition(:, 2));
extRectYmax = max(obj.WFSposition(:, 2));
octagonRectPos = [extRectXmin, extRectYmin, extRectXmax - extRectXmin, extRectYmax - extRectYmin];
gridXLength = octagonRectPos(3)*marginRatio;
gridYLength = octagonRectPos(4)*marginRatio;
centerX = (extRectXmax + extRectXmin)/2;
centerY = (extRectYmax + extRectYmin)/2;
gridMinX = centerX - gridXLength/2;
gridMaxX = centerX + gridXLength/2;
gridMinY = centerY - gridYLength/2;
gridMaxY = centerY + gridYLength/2;
xLim = [gridMinX, gridMaxX ]; yLim = [gridMinY, gridMaxY];
x = linspace(xLim(1), xLim(2), numPointsX);
y = linspace(yLim(1), yLim(2), numPoinstY);
z = 0;
[X, Y, Z] = ndgrid(x, y, z);
recPositions = [X(:), Y(:), Z(:)];
% Comment later
recPositions = [centerX, centerY, 0];

% Positions of the noise source
% Quarter of a circle
numPointsPerQuarter = 4;
radius = [5 10];
numCircles = numel(radius);
alpha = linspace(0, pi/2, numPointsPerQuarter)';
xOctagon = obj.WFSposition(:, 1);
yOctagon = obj.WFSposition(:, 2);
centreX = (max(xOctagon) + min(xOctagon))/2;
centreY = (max(yOctagon) + min(yOctagon))/2;
x = centreX + repmat(radius, numPointsPerQuarter, 1).*repmat(cos(alpha), 1, numCircles);
y = centreY + repmat(radius, numPointsPerQuarter, 1).*repmat(sin(alpha), 1, numCircles);
NSpositions = [x(:), y(:), zeros(numel(x), 1)];
% Comment later
NSpositions = [centreX + 5, centreY, 0];

% Frequencies
freqs = [440]; numFreqs = length(freqs);

% Room characteristics and impulse response of chamber
numReverbTime = 2;
beta = linspace(0, 1, numReverbTime); % Average reflection coefficient of the walls of the chamber
WFS_AcPath_previously_calculated = false;
NS_AcPath_previously_calculated = true;
appendFreeSpaceAcPaths = false;
% Comment later
beta = 0;

% WFS options
frequencyCorrection = true;

% Simulation options
timeDomainActive = true;
fakeTimeProcessing = false;
frequencyDomainActive = true;
automaticLengthModification = false;
predefSignals = true;
saveSignals = true;

durSign = 1; % Duration of tone for time processing
t = (0:ceil(durSign*obj.Fs)-1)/obj.Fs;
NSsignal = chirp(t, 20, durSign, 940);

SetupParametersScript
AcousticPathCalculationScript
simulationScript

numSamp = size(recNS_signals, 2);
t = (0:numSamp-1)/fs;

ax = axes(figure);
plot(ax, t, recNS_signals(1,:))
% axFC.XLabel.String = 'Time (s)';
% axFC.YLabel.String = 'Signal (arbitrary units)';
% axFC.Title.String = 'With Frequency Filter';

%%

obj = SimulationController;

% Constants
c = 340; % Sound velocity (m/s)
fs = 44100; % Sample frequency (samples/s)
WFSarrayOffset = [0.5, 1.5, 1.5]; % [x, y, z] coordinates. Useful for generating acoustic path IR.

% Noise source coefficient
amplitude = 1;
phase = 0;

% Filter variables for the time WFS filter.
WFSfilterLength = 22050;
% Creation of frequency filters with different orders.
magnFiltOrder = 2.^(12);
hilbertFiltOrder = 2.^(12);
numFreqFilters = length(magnFiltOrder);

freqFilters = cell(numFreqFilters, 1);
freqFiltDelays = zeros(numFreqFilters, 1);
for k = 1:numFreqFilters
    [freqFilter, delay] = getFrequencyFilter( magnFiltOrder(k), hilbertFiltOrder(k), fs );    
    freqFilters{k} = freqFilter;
    freqFiltDelays(k) = delay;
end

extRectXmin = min(obj.WFSposition(:, 1));
extRectXmax = max(obj.WFSposition(:, 1));
extRectYmin = min(obj.WFSposition(:, 2));
extRectYmax = max(obj.WFSposition(:, 2));
centerX = (extRectXmax + extRectXmin)/2;
centerY = (extRectYmax + extRectYmin)/2;
recPositions = [centerX, centerY, 0; centerX, centerY + 1, 0];

NSpositions = [3.35 -0.2 0];

% Frequencies
freqs = [440]; numFreqs = length(freqs);

% Room characteristics and impulse response of chamber
WFS_AcPath_previously_calculated = false;
NS_AcPath_previously_calculated = true;
appendFreeSpaceAcPaths = false;
% Comment later
beta = 0;

% WFS options
frequencyCorrection = true;

% Simulation options
timeDomainActive = true;
fakeTimeProcessing = false;
frequencyDomainActive = true;
automaticLengthModification = false;
predefSignals = true;
saveSignals = true;

durSign = 1; % Duration of tone for time processing
t = (0:ceil(durSign*obj.Fs)-1)/obj.Fs;
NSsignal = chirp(t, 20, durSign, 940);

SetupParametersScript
AcousticPathCalculationScript
simulationScript

numSamp = size(recNS_signals, 2);
t = (0:numSamp-1)/fs;

ax = axes(figure);
plot(ax, t, recNS_signals(1,:))
