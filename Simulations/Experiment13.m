%% Experiment 13

% Script used to generate the graphs for the section of the TFM report
% about non-free-space conditions.

%% Preamble
pathSetUp;

imagesPath = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\TFM\Img\';

dataPathName = [globalPath, 'Data\'];
ID = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

%% Simulation of different reflection coefficient (beta) values

%%% Simulation
if ~exist('obj', 'var') || ~isvalid(obj)
    obj = SimulationController;
    obj.WFSToolObj.fig.HandleVisibility = 'off';
end

% Constants
WFSfilterLength = 22050;
zPos = 1.65;
WFSarrayOffset = [0.46 2.21 zPos]; % [x, y, z] coordinates. Useful for generating acoustic path IR.
roomDim = [4.48, 9.13, 2.64];
fs = 44100;
c = 340;

% Frequency filters
magnFiltOrder = 2^10;
hilbertFiltOrder = 2^13;
[freqFilter, freqFiltDelay] = getFrequencyFilter(magnFiltOrder, hilbertFiltOrder, fs);
freqFilters = {freqFilter};
freqFiltDelays = freqFiltDelay;

% Noise source coefficient
amplitude = 1;
phase = 0;

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

% Frequencies
freqs = 0:10:1000;

% % Signal
% freqs = 0;
% durSign = 1; % Duration of tone for time processing
% t = (0:ceil(durSign*fs)-1)/fs;
% NSsignal = chirp(t, 20, durSign, 940);
% predefSignals = true;
% saveSignals = true;

% Room characteristics and impulse response of chamber
beta = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.5]; % Average reflection coefficient of the walls of the chamber
WFS_AcPath_previously_calculated = true;
NS_AcPath_previously_calculated = true;
appendFreeSpaceAcPaths = false;

% WFS options
frequencyCorrection = true;
attenuationType = 'Ruben';

% Simulation options
timeDomainActive = true;
fakeTimeProcessing = true;
frequencyDomainActive = false;
automaticLengthModification = false;

SetupParametersScript
AcousticPathCalculationScript
simulationScript

%%% Analysis
% recNS_signalsAux = pointWiseExtend( recNS_signals, rec_signals );
% recWFS_signals = rec_signals - recNS_signalsAux;
% 
% f = 0:1000;
% oper = @(x) freqz(x, 1, f, fs);
% recNS = oneDimOperOverMultiDimArray( oper, recNS_signals, 2);
% recWFS = oneDimOperOverMultiDimArray( oper, recWFS_signals, 2);

% % In order to calculate the global cancellation, it is convenient to create
% % a structure array
% ss = repmat(s(1), [numel(fsel), numNSpos]);
% for fIter = 1:numel(fsel)
%     for ns = 1:numNSpos
%         ss(fIter, ns).recCoef = recNS(:, fIter, 1, ns) + recWFS(:, fIter, 1, ns);
%         ss(fIter, ns).recNScoef = recNS(:, fIter, 1, ns); 
%         ss(fIter, ns).recWFScoef = recWFS(:, fIter, 1, ns);
%         ss(fIter, ns).Frequency = fsel(fIter);
%     end
% end

[sExt, corrFactInd, corrFactGlob, attenInd, attenGlob] =...
    SimulationController.addCancellationParametersToStructure(s);

freqEdges = 0:20:1000;
attenEdges = -20:5;
vecs = {1, freqs, 1:numNSpos, beta};
axCanc = gobjects(numFreqFilters, 1);
for k = 1:numReverbTime
    [~, attenGlobCurrent] = filterArrayForRepresentation(vecs, attenGlob, [2 3], 'nonIndepDimIndices', [1, k]);
    ax = histogram2D(10*log10(attenGlobCurrent), 1, freqs, freqEdges, attenEdges);
    ax.XLabel.String = 'Frequency (Hz)';
    ax.YLabel.String = 'Attenuation (dB)';
    ax.Title.String = ['$\beta = ', num2str(beta(k)), '$'];
    ax.Title.Interpreter = 'latex';
    colorbar(ax);
    axCanc(k) = ax;
end


% % Print
% sel = 2:5;
% for k = 1:numel(sel)
%     printfig(axCanc(sel(k)).Parent, imagesPath, ['Experiment13_globalAttenReflCoef_', num2str(beta(sel(k)))], 'eps')
% end
