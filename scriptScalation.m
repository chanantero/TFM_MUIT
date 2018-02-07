%% Why the default WFS calculation needs to be scaled?
% Study with simulation

%% Preamble
pathSetUp;

dataPathName = [globalPath, 'Data\'];
ID = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

%% System set up.
obj = SimulationController;

obj.NSposition = [3.35 -0.2 0]; % Assumed real position
obj.amplitude = 0.5;
obj.phase = 0;
obj.frequency = 440;

% Microphone positions
% Rectangular grid
extRectXmin = min(obj.WFSposition(:, 1));
extRectXmax = max(obj.WFSposition(:, 1));
extRectYmin = min(obj.WFSposition(:, 2));
extRectYmax = max(obj.WFSposition(:, 2));
octagonRectPos = [extRectXmin, extRectYmin, extRectXmax - extRectXmin, extRectYmax - extRectYmin];
gridXLength = octagonRectPos(3)*0.6;
gridYLength = octagonRectPos(4)*0.6;
centerX = (extRectXmax + extRectXmin)/2;
centerY = (extRectYmax + extRectYmin)/2;
gridMinX = centerX - gridXLength/2;
gridMaxX = centerX + gridXLength/2;
gridMinY = centerY - gridYLength/2;
gridMaxY = centerY + gridYLength/2;
xLim = [gridMinX, gridMaxX ]; yLim = [gridMinY, gridMaxY];
x = linspace(xLim(1), xLim(2), 20);
y = linspace(yLim(1), yLim(2), 20);
z = 0;
[X, Y, Z] = ndgrid(x, y, z);
grid = [X(:), Y(:), Z(:)];
obj.microPos = grid;

% Acoustic paths
obj.setAcousticPaths('NS', 'theoretical', 'WFS', 'theoretical');

%% Grid of points for the noise source
% % Rectangular grid
% gridMinX = 3;
% gridMaxX = 4;
% gridMinY = -1;
% gridMaxY = 1;
% xLim = [gridMinX, gridMaxX ]; yLim = [gridMinY, gridMaxY];
% numPointsX = 10; numPointsY = 10;
% x = linspace(xLim(1), xLim(2), numPointsX);
% y = linspace(yLim(1), yLim(2), numPointsY);
% z = 0;
% [X, Y, Z] = ndgrid(x, y, z);
% NSpositions = [X(:), Y(:), Z(:)];

xOctagon = obj.WFSposition(1:8:end, 1);
yOctagon = obj.WFSposition(1:8:end, 2);
centreX = (max(xOctagon) - min(xOctagon))/2;
centreY = (max(yOctagon) - min(yOctagon))/2;
F = logspace(log10(1.5), log10(100), 1);
x = centreX + (xOctagon - centreX)*F;
y = centreY + (yOctagon - centreY)*F;

numOct = numel(F);
numPointsPerOct = size(xOctagon, 1);

x = x(:);
y = y(:);

NSpositions = [x, y, zeros(size(x))];
numNSpos = size(NSpositions, 1);

obj.cancelResults = [];
for p = 1:numNSpos
    fprintf('%d/%d\n', p, numNSpos);
    
    obj.NSposition = NSpositions(p, :); % Assumed real position
    obj.setAcousticPaths('NS', 'theoretical', 'WFS', 'theoretical');
    
    % a) Perform WFS calculation without optimization
    obj.cancel();
        
    % b) Optimization with default theoretical acoustic paths and AllTogether
%     obj.cancel({'Loudspeakers'}, false, {'Theoretical'}, {'AllTogether'}, false);
    
    % c) Optimization with default grid of theoretical acoustic paths and
    % no restrictions (Independent).
%     obj.cancel({'Loudspeakers'}, false, {'Theoretical'}, {'Independent'}, false);
    
    % d) Optimization with theoretical acoustic paths where the microphones
    % are and AllTogether
%     obj.cancel({'Loudspeakers'}, false, {'Theoretical'}, {'AllTogether'}, false, 'testPoints', obj.microPos);
    
    % e) Optimization with theoretical acoustic paths where the microphones
    % are and Independent
%     obj.cancel({'Loudspeakers'}, false, {'Theoretical'}, {'Independent'}, false, 'testPoints', obj.microPos);
    
    % f) Optimization with current aocustic paths and AllToghether option.
    % If we use theoretical acoustic paths as current acoustic paths, this result should be the same 
    % as in d)
    obj.cancel({'Loudspeakers'}, false, {'Current'}, {'AllTogether'}, false);
    
    % g) Optimization with no restrictions. 
    % If we use theoretical acoustic paths as current acoustic paths, this result should be the same 
    % as in e)
%     obj.cancel({'Loudspeakers'}, true, {'Current'}, {'Independent'}, false);
    
    % h) Optimization with no restrictions. Special case.
    % Use solution of case f) as initial estimation of the solution with
    % magnitude constraint.
    obj.WFScoef = obj.cancelResults(end).WFScoef;
    obj.WFSToolObj.WFS_optimisation('SourceFilter', 'Loudspeakers', 'AcousticPath', 'Current', 'Grouping', 'Independent', 'maxAbsoluteValueConstraint', true, 'zerosFixed', false);
    obj.WFSToolObj.simulate();
    simulField = obj.WFSToolObj.simulField;
    sAux = obj.generateBasicExportStructure();
    obj.cancelResults = [obj.cancelResults; sAux];
end

s = obj.cancelResults;
s = reshape(s, [3, numPointsPerOct, numOct]);

% save([globalPath, 'Data/differentOcts_', ID, '.mat'], 's');

%% Visualization: 2D map, case by case

% Extract 2D map
fMap = figure;
ax = copyobj(obj.ax, fMap);
ax.Units = 'Normalized';
ax.OuterPosition = [0.5, 0, 0.5, 1];
colormap(ax, 'jet')
colorbar(ax)
% Create histogram axes
axHist = axes(fMap, 'Units', 'normalized', 'OuterPosition', [0 0 0.5 1]);

% Format structure
for p = 1:numel(s)
    s(p).NScoef = [s(p).NSRcoef; s(p).NSVcoef];
    s(p).NSposition = [s(p).NSRposition; s(p).NSVposition];
end

% Create simulationViewer object
objVis = simulationViewer(ax, s, axHist);

% % Export 2D map
% printfig(fMap, 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\', 'prueba2DMap', 'pdf');

%% Visualization: global cancellation
Cglobal = zeros(size(s));
for k = 1:numel(s)
    % Calculate global cancellation
    Cglobal(k) = sum(abs(s(k).recCoef).^2)/sum(abs(s(k).recNScoef).^2);
end

Cg_dB = 10*log10(Cglobal);

% Visualize
numPointsPerOct = size(s, 2); numOct = size(s, 3);
visualObj = animation({1:size(s, 1), 1:numPointsPerOct, 1:numOct},...
    {Cg_dB}, {'Type', 'Points', 'Octogon'}, {'Cancellation'}, [], []);

%% Visualization: correction factor AllTogether
corrFact1 = zeros(numPointsPerOct, numOct);
for o = 1:numOct
    for p = 1:numPointsPerOct
        % Calculate scaling correction factor in case b) and d)
            % Find element that is not 0
        defaultWFScoef = s(1, p, o).WFScoef;
        ind = find(defaultWFScoef ~= 0, 1, 'first');
            % Calculate
        corrFact1(p, o) = mean(s(2, p, o).WFScoef(ind)./defaultWFScoef(ind));
    end
end

ax = axes(figure);
ax.Title.Interpreter = 'latex';
hold on
data = corrFact1;
cmap = colormap('lines');
for k = 1:size(data, 2)
xData = real(data(:, k));
yData = imag(data(:, k));
scat = scatter(ax, xData, yData, 10, cmap(k, :), 'filled');
% plot(ax, xData, yData, '-', 'Color', cmap(k, :));
end
% scat.CData = 1:numPointsPerOct;
% colorbar
% ax.XLim = [-max(abs(xData(:))), max(abs(xData(:)))];
% ax.YLim = [-max(abs(yData(:))), max(abs(yData(:)))];
maxAbs = max(abs(data(:)));
ax.XLim = [-maxAbs, maxAbs]; ax.YLim = [-maxAbs, maxAbs];
ax.DataAspectRatio = [1 1 3];

meanData = mean(data);
stdData = std(data);
normStd = stdData./abs(meanData);

meanAbs = mean(abs(data));
stdAbs = std(abs(data));
normStdAbs = stdAbs./meanAbs;

meanPhase = mean(rad2deg(angle(data)));
stdPhase = std(rad2deg(angle(data)));
normStdPhase = stdPhase./meanPhase;

str = cell(size(data, 2), 1);
for k = 1:size(data, 2)
    str{k} = sprintf('$|\\overline{x}| = %.2g$, angle$(\\overline{x}) = %.0f$, $\\sigma_{n} = %.2g$, $\\sigma_{n, abs} = %.2g$, $\\sigma_{n, phase} = %.2g$', meanAbs(k), meanPhase(k), normStd(k), normStdAbs(k), normStdPhase(k));
end
l = legend(ax, str);
l.Interpreter = 'latex';

ax = axes(figure);
plot(ax, abs(data))
plot(ax, rad2deg(unwrap(angle(data))))
