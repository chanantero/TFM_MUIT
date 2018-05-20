%% Experiment 8

% Similar to Experiment7.m, but using mainly time processing, even if it's
% slower than the frequency processing.

% The analysis won't be so exhaustive, since the purpose is just proving
% that the time processing works and that the analysis processing would be
% a valid option for analysis because the correspondence between the two
% has been proven.

%% Preamble
pathSetUp;

imagesPath = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\TFM\Img\';

dataPathName = [globalPath, 'Data\'];
ID = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

%% System set up.
% obj = SimulationController;

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
% magnFiltOrder = 2.^(12);
% hilbertFiltOrder = 2.^(12);
% numFreqFilters = length(magnFiltOrder);
% 
% freqFilters = cell(numFreqFilters, 1);
% freqFiltDelays = zeros(numFreqFilters, 1);
% for k = 1:numFreqFilters
%     [freqFilter, delay] = getFrequencyFilter( magnFiltOrder(k), hilbertFiltOrder(k), fs );    
%     freqFilters{k} = freqFilter;
%     freqFiltDelays(k) = delay;
% end

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
NSpositions = [centreX + 5, centerY, 0];

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

% Duration of tone for time processing
durSign = 1;

%% Setup first parameters
SetupParametersScript

%% Pre-calculate impulse responses
% Since it would take a lot of time to calculate the IR during the
% simulations, we calculate it previously and save it to a .mat. 

AcousticPathCalculationScript
% WFSposition = obj.WFSposition;
% save([dataPathName, 'AcPathResponsesWFS_', ID, '.mat'], 'r', 'roomDim', 'numSamp', 'WFSposition', 'freqs', 'WFS_IR', 'WFS_FR')
% save([dataPathName, 'AcPathResponsesNS_', ID, '.mat'], 'r', 'roomDim', 'numSamp', 'NSpositions', 'NS_IR', 'freqs', 'NS_FR')

%% Simulation
simulationScript

%% Analysis of results
N = size(rec_signal, 2);
taux = (0:N-1)/fs;

ax = axes(figure);
plot(ax, taux, x, taux, recNS_signal(1,:)', taux, rec_signal(1,:)')
