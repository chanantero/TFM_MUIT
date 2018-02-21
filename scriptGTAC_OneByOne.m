%% GTAC anechoic chamber
% Description.
% In this script the official acoustic paths of the anechoic chamber are
% used to simulate different scenarios. The code for reproducing some
% concrete cases is used.

%% Preamble
pathSetUp;

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
% [acousticPath, microphonePositions] = ImportImpulseResponseGTAC(frequency)
load([dataPathName, 'acousticPathsGTAC_440.mat'])

% Receivers
obj.microPos = microphonePositions;

% GTAC's official acoustic paths
acPathWFSarrayStruct.acousticPaths = acousticPath;
acPathWFSarrayStruct.frequencies = 440;

obj.setAcousticPaths('NS', 'theoretical', 'WFS', acPathWFSarrayStruct);

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


%% Diferencia de fase del caso real y el caso teórico

numScen = numel(s);
numMicro = size(obj.microPos, 1);
WFS2NSratio = zeros(numMicro, numScen);
for k = 1:numScen
    WFS2NSratio(:, k) = s(k).recWFScoef./s(k).recNScoef;
end

phaseShift = angle(WFS2NSratio);
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

%% Evolución del campo en función de la distancia

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

% Magnitude
    % Specific Loudspeaker
ax = axes(figure, 'NextPlot', 'Add');
indLoud = 50;
a = abs(WFSfield(:, indLoud));
d = distances(:, indLoud);
scatter(ax, distances(:, indLoud), abs(NSfield(:, indLoud)), '.')
scatter(ax, distances(:, indLoud), abs(WFSfield(:, indLoud)), '.')

    % Global magnitude
ax = axes(figure, 'NextPlot', 'Add');
scatter(ax, distances(:), abs(NSfield(:)), '.')
scatter(ax, distances(:), abs(WFSfield(:)), '.')

    % Fitting
models = {'exp1', 'exp2', 'power1', 'power2', 'a*exp(-b*x)+c'};
numModels = numel(models);
data = repmat(struct('x', [], 'y', []), numLoudspeakers + 1, 1);
for indLoud = 1:numLoudspeakers
    data(indLoud).x = distances(:, indLoud);
    data(indLoud).y = abs(WFSfield(:, indLoud));
end
data(numLoudspeakers + 1).x = distances(:);
data(numLoudspeakers + 1).y = abs(WFSfield(:));
[params, gofs] = fitInterface( data, models );

    % Representation of fitted RMSE
ax = axes(figure, 'NextPlot', 'Add');
for m = 1:numModels
    plot(ax, [gofs(:, m).rmse])
end
legend(models)

ax = axes(figure);
plot(ax, [params{1}.b])

% Phase
    % Specific Loudspeaker
ax = axes(figure, 'NextPlot', 'Add');
indLoud = 50;
scatter(ax, distances(:, indLoud), rad2deg(angle(NSfield(:, indLoud))), '.')
scatter(ax, distances(:, indLoud), rad2deg(angle(WFSfield(:, indLoud))), '.')

    % Global phase
ax = axes(figure, 'NextPlot', 'Add');
scatter(ax, distances(:), rad2deg(angle(NSfield(:))), '.')
scatter(ax, distances(:), rad2deg(angle(WFSfield(:))), '.')

%% Optimización de la posición de la malla
% Quizás la diferencia entre el caso ideal y los casos reales se debe a un
% mal posicionamiento de la malla. Podemos optimizar la posición para
% comprobar qué offset hay que añadirle para que el parecido entre el caso
% teórico y el caso real (con los caminos acústicos oficiales) se parezcan
% lo máximo.
% ¿Cómo mido el parecido? El escalado no es importante, así que siempre se
% optimizará.

% Select the loudspeaker for which you want to optimize
optOffsets = zeros(obj.numWFS, 3);
n = 0;
for indLoud = 1:obj.numWFS
    fprintf(repmat('\b', 1, n));
    msg = sprintf('%d/%d', indLoud, obj.numWFS);
    fprintf(msg);
    n = numel(msg);
    
    obj.WFScoef(:) = 0;
    obj.WFScoef(indLoud) = 1;
    obj.NScoef = 1;
    obj.NSposition = obj.WFSposition(indLoud,:);

    % % Optimize. It didn't work
    [NSposOpt, fVal] = obj.findBestNoiseSourcePosition();
    NSOffsetOpt = NSposOpt - obj.WFSposition(indLoud,:);
    
    optOffsets(indLoud, :) = NSOffsetOpt;
end

ax = axes(figure);
scat = scatter3(ax, optOffsets(:, 1), optOffsets(:, 2), optOffsets(:, 3), '.');
scat.CData = 1:numel(scat.XData);
colorbar;
ax.XLabel.String = 'x';
ax.YLabel.String = 'y';
ax.Title.String = 'Loudspeaker Offset';
ax.DataAspectRatio = [1 1 1];
maxValue = max(abs(optOffsets(:)));
ax.XLim = [-maxValue, maxValue];
ax.YLim = [-maxValue, maxValue];

%% Pasamos a modo manual
xOffsetLim = [-0.3 -0.1];
yOffsetLim = [0.05 0.1];
numPointsX = 10;
numPointsY = 10;
xVec = linspace(xOffsetLim(1), xOffsetLim(2), numPointsX);
yVec = linspace(yOffsetLim(1), yOffsetLim(2), numPointsY);
[x, y] = ndgrid(xVec, yVec);
x = x(:);
y = y(:);
gridOffset = [x, y, zeros(numel(x), 1)];
numPointsGrid = size(gridOffset, 1);
gridPoints = repmat(obj.WFSposition(indLoud,:), numPointsGrid, 1) + gridOffset;

obj.NScoef = 1;
obj.cancelResults = [];
n = 0;
for p = 1:numPointsGrid
    fprintf(repmat('\b', 1, n));
    msg = sprintf('%d/%d', p, numPointsGrid);
    fprintf(msg);
    n = numel(msg);
    
% Place the noise source at the same position as the loudspeaker
obj.NSposition = gridPoints(p, :);

% Update the acoustic path for the noise source so it has a theoric
% acoustic path
obj.setAcousticPaths('NS', 'theoretical');

% Simulate the produced field in the receiver positions
obj.WFSToolObj.prepareOptimization();
obj.WFScoef(indLoud) = 1;
obj.WFSToolObj.WFS_optimisation('SourceFilter', 'NoFilter',...
    'AcousticPath', 'Current', 'Grouping', 'AllTogether');
obj.WFSToolObj.prepareSimulation();
obj.WFSToolObj.simulate();

% Save the information of this simulation
sAux = obj.generateBasicExportStructure();
obj.cancelResults = [obj.cancelResults; sAux];
end
fprintf('\n')

s = obj.cancelResults;

%% Visualization: global cancellation
Cglobal = zeros(size(s));
for k = 1:numel(s)
    % Calculate global cancellation
    Cglobal(k) = sum(abs(s(k).recCoef).^2)/sum(abs(s(k).recNScoef).^2);
end

Cg_dB = 10*log10(Cglobal);

% Visualize
simulFieldFormatted = mergeAndPermute(Cglobal, {3, [1 2]}, true, [numPointsX, numPointsY]);
visualObj = animation({xVec, yVec},...
    {simulFieldFormatted}, {'x', 'y'}, {'Power'}, [], []);

%% Non-conclusive analysis

% Visualize results
microPos = obj.microPos;
numMicro = size(microPos, 1);

distances = zeros(numMicro, numPointsGrid);
WFSfield = zeros(numMicro, numPointsGrid);
NSfield = zeros(numMicro, numPointsGrid);
field = zeros(numMicro, numPointsGrid);
for p = 1:numPointsGrid
    distances(:, p) = sqrt(sum((microPos - repmat(s(p).NSRposition, numMicro, 1)).^2, 2));
    WFSfield(:, p) = s(p).recWFScoef;
    NSfield(:, p) = s(p).recNScoef;
    field(:, p) = s(p).recCoef;
end
distances0 = sqrt(sum((microPos - repmat(obj.WFSposition(indLoud,:), numMicro, 1)).^2, 2));


% Magnitude
ax = axes(figure);
scat = scatter(ax, distances(:, p), abs(WFSfield(:, 1)), '.');
ax.XLim = [0 5];
for p = 1:numPointsGrid
    ax.Title.String = ['[', num2str(s(p).NSRposition), ']'];
    scat.XData = distances(:, p);
    pause(0.2)
end

% Phase
ax = axes(figure, 'NextPlot', 'Add');
ax.XLim = [0 5];
scat = scatter(ax, distances(:, 1), rad2deg(angle(WFSfield(:, 1))), '.');
scatNS = scatter(ax, distances(:, 1), rad2deg(angle(NSfield(:, 1))), '.');
for p = 1:numPointsGrid
    scat.XData = distances(:, p);
    scatNS.XData = distances(:, p);
    scatNS.YData = rad2deg(angle(NSfield(:, p)));
    ax.Title.String = ['[', num2str(s(p).NSRposition), ']'];
    
    pause(0.1)
end

ax = axes(figure, 'NextPlot', 'Add');
scat = scatter(ax, 1:obj.numMicro, rad2deg(angle(NSfield(:, 1))), '.');
for p = 1:numPointsGrid
    scat.YData = rad2deg(angle(WFSfield(:, 1)./NSfield(:, p)));
    ax.Title.String = ['[', num2str(s(p).NSRposition), ']'];
    pause(0.1)
end

phases = angle(WFSfield(:, 1));
[~, ind] = sort(distances0);
orderedPhases = phases(ind);
orderedPhases = unwrap(orderedPhases);
ax = axes(figure);
plot(ax, distances0(ind), orderedPhases)

