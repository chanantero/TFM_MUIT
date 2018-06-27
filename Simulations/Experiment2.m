%% Experiment 2. GTAC official acoustic paths.
% The original code is taken from scriptGTAC.m

%% Preamble
pathSetUp;

imagesPath = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\';

dataPathName = [globalPath, 'Data\'];
ID = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

%% System set up.
obj = SimulationController;

obj.NSposition = [3.35 -0.2 0]; % Assumed real position
obj.amplitude = 1;
obj.phase = 0;
obj.frequency = 440;

% Load information about receiver positions and acoustic paths. The
% information was generated with the command:
% [acousticPath, microphonePositions] = importImpulseResponseGTAC(frequency)
load([dataPathName, 'acousticPathsGTAC_440.mat'])

% Receivers
obj.microPos = microphonePositions;

% GTAC's official acoustic paths
acPathWFSarrayStruct.acousticPaths = acousticPath;
acPathWFSarrayStruct.frequencies = 440;

obj.setAcousticPaths('NS', 'theoretical', 'WFS', acPathWFSarrayStruct);

%% Cancellation

% Circles
numPointsPerCircle = 60;
radius = [5 50 5000];
numCircles = numel(radius);
alpha = linspace(0, 2*pi, numPointsPerCircle + 1); alpha = alpha(1:end-1)';
xOctagon = obj.WFSposition(:, 1);
yOctagon = obj.WFSposition(:, 2);
centreX = (max(xOctagon) + min(xOctagon))/2;
centreY = (max(yOctagon) + min(yOctagon))/2;
x = centreX + repmat(radius, numPointsPerCircle, 1).*repmat(cos(alpha), 1, numCircles);
y = centreY + repmat(radius, numPointsPerCircle, 1).*repmat(sin(alpha), 1, numCircles);
yOctagon = obj.WFSposition(1:8:end, 2);

x = x(:);
y = y(:);

NSpositions = [x, y, zeros(size(x))];
numNSpos = size(NSpositions, 1);

obj.cancelResults = [];
for p = 1:numNSpos
    fprintf('%d/%d\n', p, numNSpos);
    
    obj.NSposition = NSpositions(p, :); % Assumed real position
    obj.setAcousticPaths('NS', 'theoretical');
    
    obj.cancel();   
    obj.cancel({'NoFilter'}, false, {'Current'}, {'AllTogether'}, false);
    obj.cancel({'NoFilter'}, false, {'Current'}, {'Independent'}, false);
end

s = obj.cancelResults;
s = reshape(s, [3, numPointsPerCircle, numCircles]);

%% Visualize global cancellation

Cglobal = zeros(size(s));
for k = 1:numel(s)
    % Calculate global cancellation
    Cglobal(k) = sum(abs(s(k).recCoef).^2)/sum(abs(s(k).recNScoef).^2);
end

Cg_dB = 10*log10(Cglobal);

% Generate graph for TFM report
ax = axes(figure);
Cg_dB_scal = permute(Cg_dB(2, :, :), [2, 3, 1]);
plot(ax, rad2deg([alpha; 2*pi]), [Cg_dB_scal; Cg_dB_scal(1, :)])
ax.XLabel.String = '\alpha (º)';
ax.YLabel.String = 'Global cancellation (dB)';
ax.XTick = 0:90:360;
ax.XLim = [0, 360];
legLab = cell(numCircles, 1);
for c = 1:numCircles
    legLab{c} = num2str(radius(c));
end
l = legend(ax, legLab);
title(l, 'R')

% Print graph
fontSize_axesWidth_ratio = 0.08;
fontSize = ax.Position(3) * fontSize_axesWidth_ratio;
ax.XLabel.FontUnits = 'normalized';
ax.XLabel.FontSize = fontSize;
ax.YLabel.FontUnits = 'normalized';
ax.YLabel.FontSize = fontSize;

% printfig(ax.Parent, imagesPath, 'Experiment2_globalCancDifNSpos', 'eps');

% No constraint optimization
ax = axes(figure);
Cg_dB_scal = permute(Cg_dB(3, :, :), [2, 3, 1]);
plot(ax, rad2deg([alpha; 2*pi]), [Cg_dB_scal; Cg_dB_scal(1, :)])
ax.XLabel.String = '\alpha (º)';
ax.YLabel.String = 'Global cancellation (dB)';
ax.XTick = 0:90:360;
ax.XLim = [0, 360];
legLab = cell(numCircles, 1);
for c = 1:numCircles
    legLab{c} = num2str(radius(c));
end
l = legend(ax, legLab);
title(l, 'R')

fontSize_axesWidth_ratio = 0.08;
fontSize = ax.Position(3) * fontSize_axesWidth_ratio;
ax.XLabel.FontUnits = 'normalized';
ax.XLabel.FontSize = fontSize;
ax.YLabel.FontUnits = 'normalized';
ax.YLabel.FontSize = fontSize;
% printfig(ax.Parent, imagesPath, 'Experiment2_globalCancDifNSposNoConst', 'eps');

%% One by one loudspeaker excitation
% Visualize the response of each loudspeaker of the array and compare it's
% response to one of an ideal monopole

obj.cancelResults = [];
n = 0;
for indLoud = 1:obj.numWFSsources
    fprintf(repmat('\b', 1, n));
    msg = sprintf('%d/%d', indLoud, obj.numWFSsources);
    fprintf(msg);
    n = numel(msg);
    
% Place the noise source at the same position as the loudspeaker
obj.NSposition = obj.WFSposition(indLoud, :);

% Update the acoustic path for the noise source so it has a theoric
% acoustic path
obj.setAcousticPaths('NS', 'theoretical');

% All WFS loudspeakers should be silent except for the one selected
obj.WFScoef(:) = 0;
obj.WFScoef(indLoud) = 1;

% Simulate the produced field in the receiver positions
obj.WFSToolObj.prepareSimulation();
obj.WFSToolObj.simulate();

% Save the information of this simulation
sAux = obj.generateBasicExportStructure();
obj.cancelResults = [obj.cancelResults; sAux];
end
fprintf('\n')

s = obj.cancelResults;


%% Visualize 2D map

% Format structure
for p = 1:numel(s)
    s(p).NScoef = [s(p).NSRcoef; s(p).NSVcoef];
    s(p).NSposition = [s(p).NSRposition; s(p).NSVposition];
end

% Create simulationViewer object
objVis = simulationViewer(obj.ax, s);

% % Export 2D map
% printfig(objVis.fig, 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\', 'prueba2DMap', 'pdf');

%% Images individual loudspeaker
% Play with the visualization until you are satisfied.
indLoud = 20;
cmap = colormap(objVis.ax2Dmap); % basicColormap = [0 1 0; 1 1 1; 1 1 0]; cmap = extendColormap( basicColormap, 128 );
scat = findobj(objVis.ax2Dmap.Children, 'Tag', 'receiver');
objVis.scenInd = indLoud;

objVis.magnitude = 'NS';
objVis.representationType = 'abs';
objVis.ax2Dmap.CLim = [0 3];

cdata = scat.CData; % cdata = abs(s(indLoud).recNScoef);
indices = scaled2indexedColors(size(cmap, 1), objVis.ax2Dmap.CLim, cdata);
rgb = cmap(indices, :);
colorsIdealAbs = cellstr(rgb2hex(rgb));

objVis.representationType = 'phase';
cdata = scat.CData; % cdata = rad2deg(angle(s(indLoud).recNScoef));
indices = scaled2indexedColors(size(cmap, 1), objVis.ax2Dmap.CLim, cdata);
rgb = cmap(indices, :);
colorsIdealPhase = cellstr(rgb2hex(rgb));

objVis.magnitude = 'WFS';
objVis.representationType = 'abs';
objVis.ax2Dmap.CLim = [0 0.1];
cdata = scat.CData;
indices = scaled2indexedColors(size(cmap, 1), objVis.ax2Dmap.CLim, cdata);
rgb = cmap(indices, :);
colorsRealAbs = cellstr(rgb2hex(rgb));

objVis.representationType = 'phase';
cdata = scat.CData; % cdata = rad2deg(angle(s(indLoud).recNScoef));
indices = scaled2indexedColors(size(cmap, 1), objVis.ax2Dmap.CLim, cdata);
rgb = cmap(indices, :);
colorsRealPhase = cellstr(rgb2hex(rgb));

colors = {colorsIdealAbs, colorsIdealPhase, colorsRealAbs, colorsRealPhase};
names = {'Experiment2_loud20_IdealAbsMap', 'Experiment2_loud20_IdealPhaseMap',...
    'Experiment2_loud20_RealAbsMap', 'Experiment2_loud20_RealPhaseMap'};

extRectXmin = min(obj.WFSposition(:, 1));
extRectXmax = max(obj.WFSposition(:, 1));
extRectYmin = min(obj.WFSposition(:, 2));
extRectYmax = max(obj.WFSposition(:, 2));
octagonRectPos = [extRectXmin, extRectYmin, extRectXmax - extRectXmin, extRectYmax - extRectYmin];
marginX = 0.25; marginY = 0.25;
viewBox = [octagonRectPos(1) - marginX, octagonRectPos(2) - marginY, octagonRectPos(3) + 2*marginX, octagonRectPos(4) + 2*marginY];
WFScolor = repmat({'#000000'}, obj.numWFS, 1);
WFScolor{indLoud} = '#AA0000';
objSVG = SVGdrawer('viewBox', viewBox, 'NSpositions', double.empty(0,2),...
    'microSymbol', 'dot', 'microSize', 0.1, 'microPos', obj.microPos, 'WFScolor', WFScolor);

currentFolder = pwd;
cd(imagesPath); % Needed for inkscape to link svg files properly
for k = 1:numel(names)
    name = names{k};
    objSVG.microColor = colors{k};
    fileName = [imagesPath, name, '.svg'];
    objSVG.drawSVG(fileName);
    system(['inkscape -z "', imagesPath, name, '.svg" --export-pdf="', imagesPath, name, '.pdf"'])
end
cd(currentFolder)

%% Evolution of field as a function of distance to loudspeakaer

numLoudspeakers = numel(s);
microPos = obj.microPos;
numMicro = size(microPos, 1);

distances = zeros(numMicro, numLoudspeakers);
WFSfield = zeros(numMicro, numLoudspeakers);
NSfield = zeros(numMicro, numLoudspeakers);
for k = 1:numLoudspeakers
    distances(:, k) = sqrt(sum((microPos - repmat(obj.WFSposition(k, :), numMicro, 1)).^2, 2));
    WFSfield(:, k) = s(k).recWFScoef;
    NSfield(:, k) = s(k).recNScoef;
end

% Fitting. Useful for later.
models = {'a/x'}; %{'exp1', 'exp2', 'power1', 'power2', 'a*exp(-b*x)+c'};
data = repmat(struct('x', [], 'y', []), numLoudspeakers + 1, 1);
for indLoud = 1:numLoudspeakers
    data(indLoud).x = distances(:, indLoud);
    data(indLoud).y = abs(WFSfield(:, indLoud));
end
data(numLoudspeakers + 1).x = distances(:);
data(numLoudspeakers + 1).y = abs(WFSfield(:));
[params, gofs] = fitInterface( data, models );

% Magnitude
    % Specific Loudspeaker
ax = axes(figure, 'NextPlot', 'Add');
indLoud = 20;
a = abs(WFSfield(:, indLoud));
d = distances(:, indLoud);
dVec = 0:0.05:max(d);
% scatter(ax, distances(:, indLoud), params{1}(indLoud).a*abs(NSfield(:, indLoud)), '.')
plot(ax, dVec, params{1}(indLoud).a./dVec)
scatter(ax, distances(:, indLoud), abs(WFSfield(:, indLoud)), '.')
ax.YLim = [0 max(abs(WFSfield(:, indLoud)))];
ax.XLabel.String = 'Distance (m)';
ax.YLabel.String = 'Field';
legend(ax, {'Theoretical', 'Measured'});
printfig(ax.Parent, imagesPath, 'Experiment2_loud20_AmpByDist', 'eps');

    % Global magnitude
ax = axes(figure, 'NextPlot', 'Add');
% scatter(ax, distances(:), abs(NSfield(:)), '.')
scatter(ax, distances(:), abs(WFSfield(:)), '.')

% Phase
phaseReal = angle(WFSfield);
phaseTheo = -2*pi*obj.frequency/obj.WFSToolObj.c*distances; % It is the same as: phaseTheo = angle(NSfield);

    % Specific Loudspeaker
indLoud = 29;
phaseRealLoud = phaseReal(:, indLoud);
phaseTheoLoud = phaseTheo(:, indLoud);
distLoud = distances(:, indLoud);
% [distSort, indSort] = sort(dist);
% idealPhase = -2*pi*obj.frequency/obj.WFSToolObj.c*distSort;
% ax = axes(figure, 'NextPlot', 'Add');
% realPhase = idealPhase + wrapToPi(phase(indSort) - idealPhase);
% plot(ax, distSort, rad2deg(idealPhase))
% scatter(ax, distSort, rad2deg(realPhase), '.')
% ax.XLabel.String = 'Distance (m)';
% ax.YLabel.String = 'Phase (º)';
% legend(ax, {'Theoretical', 'Measured'});

    % Phase shift
ax = axes(figure);
scatter(ax, distLoud, rad2deg(wrapToPi(phaseLoud - phaseTheoLoud)), '.');
ax.YLim = [-180, 180];
ax.YTick = -180:60:180;
ax.YGrid = 'on';
ax.XLabel.String = 'Distance (m)';
ax.YLabel.String = 'Phase shift (º)';
% printfig(ax.Parent, imagesPath, 'Experiment2_loud20_PhaseByDist', 'eps');

    % Global phase

    % Global Phase shift
ax = axes(figure);
scatter(ax, distances(:), rad2deg(wrapToPi(phaseReal(:) - phaseTheo(:))), '.');
ax.YLim = [-180, 180];
ax.YTick = -180:60:180;
ax.YGrid = 'on';
ax.XLabel.String = 'Distance (m)';
ax.YLabel.String = 'Phase shift (º)';

     % % Other way
% ax = axes(figure, 'NextPlot', 'Add');
% phase = angle(WFSfield(:));
% dist = distances(:);
% [distSort, indSort] = sort(dist);
% idealPhase = -2*pi*obj.frequency/obj.WFSToolObj.c*distSort;
% realPhase = idealPhase + wrapToPi(phase(indSort) - idealPhase);
% plot(ax, distSort, rad2deg(idealPhase))
% scatter(ax, distSort, rad2deg(realPhase), '.')
% ax.XLabel.String = 'Distance (m)';
% ax.YLabel.String = 'Phase (º)';
% legend(ax, {'Theoretical', 'Measured'});


% Phase shift
phaseShift = wrapToPi(phaseReal - phaseTheo); % The same as: phaseShift = angle(WFSfield./NSfield);

[meanPhaseShift, ~, stdPhaseShift] = circularDistributionParameters(phaseShift);
meanPhaseShift = rad2deg(meanPhaseShift);
stdPhaseShift = rad2deg(stdPhaseShift);

ax = axes(figure);
plot(ax, 1:numScen, meanPhaseShift, 1:numScen, meanPhaseShift + stdPhaseShift, '--', 1:numScen, meanPhaseShift - stdPhaseShift, '--')
ax.XLabel.String = 'Loudspeakers';
ax.YLabel.String = 'Phase Shift (º)';

globalMean = mean(phaseShift(:));
globalStd = std(phaseShift(:));

axHist = axes(figure);
histogram(axHist, phaseShift)

%% Optimization without restrictions. What is the best cancellation we can achieve?


