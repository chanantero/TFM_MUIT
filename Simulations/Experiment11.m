%% Experiment 11.

% Based on Experiment7_rebuilding.m.
% Created on 27th June 2018, it aims at generating the graphs for the
% introduction to simulations chapter in the thesis report.

%% Preamble
pathSetUp;

imagesPath = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\TFM\Img\';

dataPathName = [globalPath, 'Data\'];
ID = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

%% System set up.
if ~exist('obj', 'var') || ~isvalid(obj)
    obj = SimulationController;
end

% Constants
c = 340; % Sound velocity (m/s)
fs = 44100; % Sample frequency (samples/s)
WFSfilterLength = 22050;
zPos = 1.65;
WFSarrayOffset = [0.46 2.21 zPos]; % [x, y, z] coordinates. Useful for generating acoustic path IR.
roomDim = [4.48, 9.13, 2.64];

% Noise source coefficient
amplitude = 1;
phase = 0;

% Filter variables for the time WFS filter. Creation of frequency filters
% with different orders.
magnFiltOrder = 2.^(12);
hilbertFiltOrder = 2.^(12);
[freqFilter, freqFiltDelay] = getFrequencyFilter( magnFiltOrder, hilbertFiltOrder, fs );    
freqFilters = {freqFilter};
freqFiltDelays = freqFiltDelay;

% Microphone positions
% Rectangular grid
marginRatio = 0.6;
numPointsX = 5;
numPoinstY = 5;
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

% Positions of the noise source
% Quarter of a circle
NSpositions = [centerX + 5, centerY, 0];

% Frequencies
freqs = 440;

% Room characteristics and impulse response of chamber
beta = 0; % Average reflection coefficient of the walls of the chamber
WFS_AcPath_previously_calculated = true;
NS_AcPath_previously_calculated = true;
appendFreeSpaceAcPaths = false;

% WFS options
frequencyCorrection = false;
attenuationType = 'Ruben';

% Simulation options
timeDomainActive = true;
fakeTimeProcessing = false;
frequencyDomainActive = true;
automaticLengthModification = false;
saveSignals = true;

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

%% Graphs for the report

% Time domain
numSamples = size(rec_signals, 2);
t = (0:numSamples - 1)/fs;
ax = axes(figure);
plot(ax, t, recNS_signals(13, :), t, rec_signals(13, :))
ax.XLim = [0.2, 0.21];
ax.XLabel.String = 'Time (s)';
% printfig(ax.Parent, imagesPath, 'Experiment11_exampleNoFreqCorr', 'eps');

[s, corrFactInd, corrFactGlob, attenInd, attenGlob] = SimulationController.addCancellationParametersToStructure(s);
attenInd(:, 1, 1, :, 13)
cancInd = -20*log10(attenInd(2, 1, 13));
corrFact = corrFactInd(2, 1, 13);

%% SVG scenario of just one receiver measure point
viewBox = [-WFSarrayOffset(1) -WFSarrayOffset(2) roomDim(1) roomDim(2)];
NSangles = atan2d(centerY - NSpositions(:,2), centerX - NSpositions(:,1));

objSVG = SVGdrawer('viewBox', viewBox, 'NSpositions', NSpositions,...
    'NSangles', NSangles, 'microSymbol', 'dot', 'microSize', 0.05,...
    'microPositions', [centerX, centerY]);

name = 'Experiment11_scheme_oneReceiver';
objSVG.drawSVG([imagesPath, name, '.svg']);

currentFolder = pwd;
cd(imagesPath); % Needed for inkscape to link svg files properly
system(['inkscape -z "', imagesPath, name, '.svg" --export-pdf="', imagesPath, name, '.pdf"'])
cd(currentFolder)

%% SVG scenario of a grid of points of measure
viewBox = [-WFSarrayOffset(1) -WFSarrayOffset(2) roomDim(1) roomDim(2)];
NSangles = atan2d(centerY - NSpositions(:,2), centerX - NSpositions(:,1));

objSVG = SVGdrawer('viewBox', viewBox, 'NSpositions', NSpositions,...
    'NSangles', NSangles, 'microSymbol', 'dot', 'microSize', 0.05,...
    'microPositions', recPositions);

name = 'Experiment11_scheme_multipleReceiver';
objSVG.drawSVG([imagesPath, name, '.svg']);

currentFolder = pwd;
cd(imagesPath); % Needed for inkscape to link svg files properly
system(['inkscape -z "', imagesPath, name, '.svg" --export-pdf="', imagesPath, name, '.pdf"'])
cd(currentFolder)

%% What happens for the grid of points?

objVis = simulationViewer(obj.ax, s);
objVis.magnitude = 'Cancellation';
objVis.representationType = 'dB';
objVis.scenInd = 2;

fig = figure;
newax = copyobj(objVis.ax2Dmap, fig);
colorbar(newax)
colormap(newax, flipud(jet))
newax.Children(2).CData = -newax.Children(2).CData;
newax.CLim = [0, 4];

% printfig(fig, imagesPath, 'Experiment11_multipleReceiverCancel2Dmap', 'eps')

cancGlob = -20*log10(attenGlob(2));

%% Different positions of noise source

numPointsPerArc = 4;
radius = [3.6 4 4.4 4.8];
numArcs = numel(radius);
xOctagon = obj.WFSposition(:, 1);
yOctagon = obj.WFSposition(:, 2);
centreX = (max(xOctagon) + min(xOctagon))/2;
centreY = (max(yOctagon) + min(yOctagon))/2;
alphaMax = pi/2;
alphaMin = 0;
alpha = linspace(alphaMin, alphaMax, numPointsPerArc)';
x = centreX + repmat(radius, numPointsPerArc, 1).*repmat(cos(alpha), 1, numArcs);
y = centreY + repmat(radius, numPointsPerArc, 1).*repmat(sin(alpha), 1, numArcs);
NSpositions = [x(:), y(:), zeros(numel(x), 1)];

SetupParametersScript

%% SVG scenario for multiple measure points and noise sources
viewBox = [-WFSarrayOffset(1), -WFSarrayOffset(2), max(NSpositions(:, 1))+WFSarrayOffset(1)+0.5, max(NSpositions(:, 2))+WFSarrayOffset(2)+0.5];
NSangles = atan2d(centerY - NSpositions(:,2), centerX - NSpositions(:,1));

objSVG = SVGdrawer('viewBox', viewBox, 'NSpositions', NSpositions,...
    'NSangles', NSangles, 'microSymbol', 'dot', 'microSize', 0.05,...
    'microPositions', recPositions);

name = 'Experiment11_scheme_multipleReceiverMultipleNS';
objSVG.drawSVG([imagesPath, name, '.svg']);

currentFolder = pwd;
cd(imagesPath); % Needed for inkscape to link svg files properly
system(['inkscape -z "', imagesPath, name, '.svg" --export-pdf="', imagesPath, name, '.pdf"'])
cd(currentFolder)