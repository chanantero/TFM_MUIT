%% Why the default WFS calculation needs to be scaled?
% Study with simulation

% Description of the problem
%
%

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

%% Cancellation attempt

% Perform WFS calculation without optimization and save the WFS
% coefficients.
obj.cancel();
WFSnoOpt = obj.WFScoef;

% Perform WFS calculation with optimization and save the WFS coefficients
% The optimization should be done with the theoretical path and alltogether
% option
% Optimization options
sourceFilter = {'Loudspeakers'}; % It makes no sense to optimize a less real scenario, so use loudspeakers filter
maxAbsValCons = false; % We always want to be realistic about real constraints
acousticPathType = {'Theoretical'}; % It makes no sense to optimize with a theoric acoustic path because it depends on the parameters of the noise source, and those parameters are actually unknown. Besides, the acousic path of the loudspeakers is only known in the places where microphones have been placed.
grouping = {'AllTogether'};
zerosFixed = false;
% Optimization
obj.cancel(sourceFilter, maxAbsValCons, acousticPathType, grouping, zerosFixed);

WFSOpt = obj.WFScoef;

rel = WFSOpt(end)./WFSnoOpt(end);

abs(rel)
rad2deg(angle(rel))

%% Visualization
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
s = obj.cancelResults;
for p = 1:numel(s)
    s(p).NScoef = [s(p).NSRcoef; s(p).NSVcoef];
    s(p).NSposition = [s(p).NSRposition; s(p).NSVposition];
end

% Create simulationViewer object
objVis = simulationViewer(ax, s, axHist);

% % Export 2D map
% printfig(fMap, 'C:\Users\Rub�n\Google Drive\Telecomunicaci�n\M�ster 2� Curso 2015-2016\TFM MUIT\Documentos\Img\', 'prueba2DMap', 'pdf');

%% Grid of points for the noise source
% Rectangular grid
gridMinX = 3;
gridMaxX = 4;
gridMinY = -1;
gridMaxY = 1;
xLim = [gridMinX, gridMaxX ]; yLim = [gridMinY, gridMaxY];
x = linspace(xLim(1), xLim(2), 10);
y = linspace(yLim(1), yLim(2), 10);
z = 0;
[X, Y, Z] = ndgrid(x, y, z);
NSpositions = [X(:), Y(:), Z(:)];
numNSpos = size(NSpositions, 1);

% Optimization options
sourceFilter = {'Loudspeakers'}; % It makes no sense to optimize a less real scenario, so use loudspeakers filter
maxAbsValCons = false; % We always want to be realistic about real constraints
acousticPathType = {'Theoretical'}; % It makes no sense to optimize with a theoric acoustic path because it depends on the parameters of the noise source, and those parameters are actually unknown. Besides, the acousic path of the loudspeakers is only known in the places where microphones have been placed.
grouping = {'AllTogether'};
zerosFixed = false;

rel = zeros(numNSpos, 1);
for p = 1:numNSpos
    obj.NSposition = NSpositions(p, :); % Assumed real position
    
    % Perform WFS calculation without optimization
    obj.cancel();
    WFSnoOpt = obj.WFScoef;
    
    % Optimization
    obj.cancel(sourceFilter, maxAbsValCons, acousticPathType, grouping, zerosFixed);
    WFSOpt = obj.WFScoef;
    
    rel(p) = WFSOpt(end)./WFSnoOpt(end);
end

rel = rel(~isnan(rel) & isfinite(rel));