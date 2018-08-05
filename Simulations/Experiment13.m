%% Experiment 13

% Script used to generate the graphs for the section of the TFM report
% about non-free-space conditions.

%% Preamble
pathSetUp;

imagesPath = 'C:\Users\Rub�n\Google Drive\Telecomunicaci�n\M�ster 2� Curso 2015-2016\TFM MUIT\Documentos\TFM\Img\';

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
% Circle arc
numPointsPerArc = 4;
radius = [3.6 4 4.4 4.8];
numArcs = numel(radius);
xOctagon = obj.WFSposition(:, 1);
yOctagon = obj.WFSposition(:, 2);
centreX = (max(xOctagon) + min(xOctagon))/2;
centreY = (max(yOctagon) + min(yOctagon))/2;
x1 = centreX + WFSarrayOffset(1);
x2 = roomDim(1) - (centreX + WFSarrayOffset(1));
alphaMax = -pi/2 - asin(x1/max(radius)) + deg2rad(1);
alphaMin = -pi/2 + asin(x2/max(radius)) - deg2rad(1);
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

%%% Analysis. For when fakeTimeProcessing = false.
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

[sExt, corrFactInd, corrFactGlob, attenInd, attenGlob, corrFactAver, gainAver] =...
    SimulationController.addCancellationParametersToStructure(s);

freqEdges = 0:20:1000;
attenEdges = -20:5;
vecs = {1, freqs, 1:numNSpos, beta};
axCanc = gobjects(numFreqFilters, 1);
for k = 1:numReverbTime
    [~, attenGlobCurrent] = filterArrayForRepresentation(vecs, attenGlob, [2 3], 'nonIndepDimIndices', [1, k]);
    ax = histogram2D(10*log10(attenGlobCurrent), 1, freqs, freqEdges, attenEdges);
    ax.XLabel.String = 'Frequency (Hz)';
    ax.YLabel.String = 'Global gain (dB)';
    ax.Title.String = ['$\beta = ', num2str(beta(k)), '$'];
    ax.Title.Interpreter = 'latex';
    colorbar(ax);
    axCanc(k) = ax;
end


% % Print
% sel = [1 3 5 7];
% for k = 1:numel(sel)
%     printfig(axCanc(sel(k)).Parent, imagesPath, ['Experiment13_globalAttenReflCoef_', num2str(round(beta(sel(k))*100))], 'eps')
% end

freqEdges = 0:20:1000;
attenEdges = -20:5;
vecs = {1, freqs, 1:numNSpos, beta};
axGainAver = gobjects(numFreqFilters, 1);
for k = 1:numReverbTime
    [~, gainAverCurrent] = filterArrayForRepresentation(vecs, gainAver, [2 3], 'nonIndepDimIndices', [1, k]);
    ax = histogram2D(10*log10(gainAverCurrent), 1, freqs, freqEdges, attenEdges);
    ax.XLabel.String = 'Frequency (Hz)';
    ax.YLabel.String = 'Average gain (dB)';
    ax.Title.String = ['$\beta = ', num2str(beta(k)), '$'];
    ax.Title.Interpreter = 'latex';
    colorbar(ax);
    axGainAver(k) = ax;
end

% Print
sel = [1 3 5 7];
for k = 1:numel(sel)
    printfig(axGainAver(sel(k)).Parent, imagesPath, ['Experiment13_gainAverReflCoef_', num2str(round(beta(sel(k))*100))], 'eps')
end

%% SVG scenario of the room
viewBox = [-WFSarrayOffset(1) -WFSarrayOffset(2) roomDim(1) roomDim(2)];
NSangles = atan2d(centerY - NSpositions(:,2), centerX - NSpositions(:,1));

objSVG = SVGdrawer('viewBox', viewBox, 'NSpositions', NSpositions,...
    'NSangles', NSangles, 'microSymbol', 'dot', 'microSize', 0.05,...
    'microPositions', recPositions);

name = 'Experiment13_scheme';
objSVG.drawSVG([imagesPath, name, '.svg']);

currentFolder = pwd;
cd(imagesPath); % Needed for inkscape to link svg files properly
system(['inkscape -z "', imagesPath, name, '.svg" --export-pdf="', imagesPath, name, '.pdf"'])
cd(currentFolder)

%% GTAC listening room impulse response analysis

% freqs = 50:100:1000;
% [acousticPath, microphonePositions] = importImpulseResponseGTAC(freqs);
% save([dataPathName, 'acousticPathsGTAC_someFreq1.mat'], 'acousticPath', 'microphonePositions', 'freqs')
load([dataPathName, 'acousticPathsGTAC_someFreq1.mat'])
[numMicro, numLoud, numFreqs] = size(acousticPath);

% Calculate distances.
% Transform everything to distance (all receivers and loudspeakers
% together). It makes no sense to select individual loudspeakers. This
% simplifies things.
scenStruct = WFSToolSimple.generateScenario(96, 'orientation', 'vertical', 'originReference', 'octagonBondingBoxCorner');
WFSarrayPosition = scenStruct.loudspeakersPosition;
distances = calcDistances(microphonePositions, WFSarrayPosition);
distances = distances(:);
acPath = mergeAndPermute(acousticPath, {[1,2], 3});

% Fitting. Useful for later.
models = {'a/x'}; %{'exp1', 'exp2', 'power1', 'power2', 'a*exp(-b*x)+c'};
data = repmat(struct('x', [], 'y', []), [numFreqs, 1]);
for f = 1:numFreqs
    data(f).x = distances;
    data(f).y = abs(acPath(:, f));
end
[params, gofs] = fitInterface( data, models );

dataAll = struct('x', repmat(distances, [numFreqs, 1]), 'y', abs(acPath(:)));
[paramAll, gofsGlob] = fitInterface( dataAll, models );

% Types of representation:
% A) Fixed frequency, all distances. Histogram/scatter. x axis: distance. y axis: magnitude.
% B) All frequencies, fixed distance (pair loudspeaker - measure point). Histogram/scatter/plot. x axis: frequency. y axis: magnitude.
% C) All frequencies, all distances. Histogram/scatter. x axis: distance. y axis: magnitude.
% D) Scatter. x axis: frequency. y axis: distance. z axis: magnitude.
% E) 2D histogram. x axis: frequency. y axis: distance. z/color axis: mean/deviation/std.

% A) Fixed frequency, all distances. Histogram/scatter. x axis: distance. y axis: magnitude.
repFreqs = freqs;
numRepFreqs = length(repFreqs);
dVec = min(distances):0.05:max(distances);
distEdges = linspace(min(distances), max(distances), 30);
magEdges = linspace(0, max(abs(acPath(:))), 50);
axsA = gobjects(length(repFreqs), 1);
for f = 1:numRepFreqs
    [~, indFreq] = min(abs(repFreqs(f) - freqs));
    repArr = abs(acPath(:, indFreq));
    ax = axes(figure, 'NextPlot', 'Add');
    scatter(ax, distances, repArr, '.')
%     histogram2D(repArr, 1, distances, distEdges, magEdges, 'axes', ax);
    plot(ax, dVec, params{1}(indFreq).a./dVec, 'r', 'LineWidth', 3)
    ax.XLim = [min(distances), max(distances)];
    ax.YLim = [0 max(abs(acPath(:)))];
    ax.XLabel.String = 'Distance (m)';
    ax.YLabel.String = '$|H|$'; ax.YLabel.Interpreter = 'latex';
    ax.Title.String = sprintf('%gHz', repFreqs(f));
    legend(ax.Children(1), {'Ideal free-space'});
    axsA(f) = ax;
end
close all

% % Print
% for f = 1:numRepFreqs
%     printfig(axsA(f).Parent, imagesPath, ['Experiment13_AmpByDist_', num2str(repFreqs(f)), 'Hz'], 'eps');
% end


% B) All frequencies, fixed distance (pair loudspeaker - measure point). Histogram/scatter/plot. x axis: distance. y axis: magnitude.
% Select index of loudspeaker and microphone and calculate index of distance
indLoud = 20;
indMicro = 130;
indDist = sub2ind([numMicro, numLoud], indMicro, indLoud);
% Select index of distance and calculate index of loudspeaker and microphone 
% indDist = 20;
% [indMicro, indLoud] = ind2sub([numMicro, numLoud], indDist);

ax = axes(figure, 'NextPlot', 'Add');
a = abs(acPath(indDist, :));
plot(ax, freqs, a);
ax.XLim = [0, max(freqs)];
ax.YLim = [0 max(abs(acPath(:)))];
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = '$|H|$'; ax.YLabel.Interpreter = 'latex';
ax.Title.String = sprintf('indLoud = %d, indMicro = %d, distance = %g', indLoud, indMicro, distances(indDist));

% C) All frequencies, all distances. Histogram/scatter. x axis: distance. y axis: magnitude.
ax = axes(figure, 'NextPlot', 'Add');
histogram2D(abs(acPath), 1, distances, distEdges, magEdges, 'axes', ax);
dVec = min(distances):0.05:max(distances);
plot(ax, dVec, paramAll{1}.a./dVec, 'r', 'LineWidth', 3)
ax.XLim = [min(distances), max(distances)];
ax.YLim = [0 max(abs(acPath(:)))];
ax.XLabel.String = 'Distance (m)';
ax.YLabel.String = '$|H|$'; ax.YLabel.Interpreter = 'latex';
legend(ax.Children(1), {'Ideal free-space'});


% D) Scatter. x axis: frequency. y axis: distance. z axis: magnitude.
ax = axes(figure);
[X, Y] = ndgrid(distances, freqs);
scatter3(ax, X(:), Y(:), abs(acPath(:)))


% E) 2D histogram. x axis: frequency. y axis: distance. z/color axis: mean/deviation/std.
distEdges = linspace(min(distances), max(distances), 30);
freqStep = freqs(2) - freqs(1);
freqEdges = [freqs, freqs(end) + freqStep] - freqStep/2;
[D, F] = ndgrid(distances, freqs);
[N, ~, ~, binX, binY] = histcounts2(D, F, distEdges, freqEdges);

as = [params{1}(:).a];
[D, A] = ndgrid(distances, as);
theoAbs = A./D;
dev = abs(acPath) - theoAbs;
rep = zeros(size(N));
for x = 1:length(distEdges) - 1
    for y = 1:length(freqEdges) - 1
%         rep(x, y) = std(abs(acPath(binX == x & binY == y)));
        rep(x, y) = sqrt(mean(dev(binX == x & binY == y).^2));
    end
end

ax = axes(figure);
bar3cRub(rep, distEdges, freqEdges, ax);

% C = zeros(length(distEdges), length(freqEdges));
% C(1:end-1, 1:end-1) = rep;
% pcolor(ax, distEdges, freqEdges, C');

% Show magnitude-distance relation for the simulated impulse responses. A
% graph for each frequency.
% Trabajar con las variables WFS_FR y WFS_IR, aunque es mejor generarlas
% aqu� con rir_generator. Har�a que el c�digo fuese m�s independiente. S�,
% haz eso.


% Generate the average gain histogram for the reflection coefficient of the most similar simulated
% impulse response

