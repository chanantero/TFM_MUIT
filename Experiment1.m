%% Experiment 1
% Most of this code is from scriptScalation.m

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

% Microphone positions
% Rectangular grid
marginRatio = 0.6;
numPointsX = 20;
numPoinstY = 20;
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
grid = [X(:), Y(:), Z(:)];
obj.microPos = grid;

% Acoustic paths
obj.setAcousticPaths('NS', 'theoretical', 'WFS', 'theoretical');

%% One position cancellation
% Choose a position for the noise source.
obj.NSposition = [-1 2.5 0];

% Apply WFS calculation with global scalation
obj.setAcousticPaths('NS', 'theoretical');
obj.cancel({'NoFilter'}, false, {'Current'}, {'AllTogether'}, false);

% Create a pixel grid for the image and set it as microphone positions
marginRatioLeft = 0.6;
marginRatioRight = 0.2;
marginRatioTop = 0.2;
marginRatioBottom = 0.2;
numPointsX = 400;
numPointsY = 400;

extRectXmin = min(obj.WFSposition(:, 1));
extRectXmax = max(obj.WFSposition(:, 1));
extRectYmin = min(obj.WFSposition(:, 2));
extRectYmax = max(obj.WFSposition(:, 2));
gridMinX = extRectXmin - gridXLength*marginRatioLeft;
gridMaxX = extRectXmax + gridXLength*marginRatioRight;
gridMinY = extRectYmin - gridYLength*marginRatioBottom;
gridMaxY = extRectYmax + gridYLength*marginRatioTop;
xLim = [gridMinX, gridMaxX]; yLim = [gridMinY, gridMaxY];
x = linspace(xLim(1), xLim(2), numPointsX);
y = linspace(yLim(1), yLim(2), numPointsY);
z = 0;
[X, Y, Z] = ndgrid(x, y, z);
imageGrid = [X(:), Y(:), Z(:)];

% Simulate
% For efficiency purposes, simulate big number of receivers at a more low
% level
field = obj.WFSToolObj.simulTheo.calculateTheoricField(imageGrid);
fieldNS = field(:, 1);
fieldWFS = field(:, 2);
field = sum(field, 2);

% Print image
viewBox = [xLim(1), yLim(1), xLim(2) - xLim(1), yLim(2) - yLim(1)];

basicColormap = [0 1 0; 1 1 1; 1 1 0];
colorMap_R = interp1(1:3, basicColormap(:,1), linspace(1, 3, 128)');
colorMap_G = interp1(1:3, basicColormap(:,2), linspace(1, 3, 128)');
colorMap_B = interp1(1:3, basicColormap(:,3), linspace(1, 3, 128)');
extendedColormap = [colorMap_R, colorMap_G, colorMap_B];

cdata = real(field);
cdata = reshape(cdata, [numPointsX, numPointsY]);
cdata = cdata.';

% Visualize image in order to check if everything is allright
ax = axes(figure);
im = image(cdata);
im.CDataMapping = 'scaled';
ax.CLim = [-1, 1];
colormap(ax, extendedColormap);

indices = scaled2indexedColors(size(extendedColormap, 1), ax.CLim, cdata);
imwrite(indices, extendedColormap, 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\Experiment1_Example_Field.png');

name = 'Experiment1_Example';
fileName = [imagesPath, name, '.svg'];
drawWFSarrayFun( fileName, viewBox, 'NSposition', obj.NSRposition(1:2), 'backgroundFileName', 'Experiment1_Example_Field.png',...
    'microPosition', obj.microPos);

currentFolder = pwd;
cd(imagesPath); % Needed for inkscape to link svg files properly
system(['inkscape -z "', imagesPath, name, '.svg" --export-pdf="', imagesPath, name, '.pdf"'])
cd(currentFolder)

%% Grid of points for the noise source
% Create image with these noise source positions

    % Grid of point: circles
numPointsPerCircle = 20;
radius = 5;
numCircles = numel(radius);
alpha = linspace(0, 2*pi, numPointsPerCircle + 1); alpha = alpha(1:end-1)';
xOctagon = obj.WFSposition(:, 1);
yOctagon = obj.WFSposition(:, 2);
centreX = (max(xOctagon) + min(xOctagon))/2;
centreY = (max(yOctagon) + min(yOctagon))/2;
x = centreX + repmat(radius, numPointsPerCircle, 1).*repmat(cos(alpha), 1, numCircles);
y = centreY + repmat(radius, numPointsPerCircle, 1).*repmat(sin(alpha), 1, numCircles);

x = x(:);
y = y(:);

NSpositions = [x, y, zeros(size(x))];

NSposRel = repmat([centreX, centreY], size(NSpositions, 1), 1) - NSpositions(:, 1:2);
NSangles = atan2d(NSposRel(:, 2), NSposRel(:, 1));

    % Create WFS array in SVG
viewBox = [min(x)-1, min(y)-1, max(x) - min(x) + 2, max(y) - min(y) + 2];
svgText = WFSarraySVG( viewBox, 'NSposition', NSpositions(:, 1:2), 'NSangle', NSangles);

    % Draw line from the center of the octagon to the circle of loudspeakers
selLoud = 2;
strLine1 = makePath('000000', 0.01, centreX, centreY, NSpositions(selLoud,1), NSpositions(selLoud,2), 'line');

    % Draw horizontal line
strLine2 = makePath('000000', 0.01, centreX, centreY, centreX + radius/2, centreY + 0, 'line');

    % Draw arc to mark the angle alpha
endAngle = NSangles(selLoud) + 180;
startAngle = 0;
strArc = makeCircumferenceArc(centreX, centreY, 1, startAngle, endAngle, '000000', 0.01);

    % Draw the symbol of the distance
dir = (NSpositions(selLoud, 1:2) - [centreX, centreY])/radius;
strDistance = makeText( '$R$', 0.1, '000000', centreX + dir(1)*radius*0.75, centreY + dir(2)*radius*0.75, 'distanceSymbol');

    % Draw the symbol for the angle
strAngle = makeText( '$\alpha$', 0.1, '000000', centreX + 1, centreY + 0.4, 'distanceSymbol');

svgText = strrep(svgText, '[Other]', [strLine1, sprintf('\n'), strLine2, sprintf('\n'),...
    strArc, sprintf('\n'), strAngle, sprintf('\n'), strDistance]);

    % Write SVG file
name = 'Experiment1_differentNSpositions';
fileName = [imagesPath, name, '.svg'];
destFile = fopen(fileName, 'w', 'n', 'UTF-8');
fwrite(destFile, svgText, 'char');
fclose(destFile);

    % Export to PDF
currentFolder = pwd;
cd(imagesPath); % Needed for inkscape to link svg files properly
system(['inkscape -z "', imagesPath, name, '.svg" --export-pdf="', imagesPath, name, '.pdf" --export-latex'])
cd(currentFolder)

%% Apply cancellation

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
end

s = obj.cancelResults;
s = reshape(s, [2, numPointsPerCircle, numCircles]);

% save([globalPath, 'Data/differentOcts_', ID, '.mat'], 's');

%% Visualization: 2D map, case by case

% Format structure
for p = 1:numel(s)
    s(p).NScoef = [s(p).NSRcoef; s(p).NSVcoef];
    s(p).NSposition = [s(p).NSRposition; s(p).NSVposition];
end

% Create simulationViewer object
objVis = simulationViewer(obj.ax, s);

% % Export 2D map
% printfig(fMap, 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\', 'prueba2DMap', 'pdf');

%% Visualization: global cancellation
Cglobal = zeros(size(s));
for k = 1:numel(s)
    % Calculate global cancellation
    Cglobal(k) = sum(abs(s(k).recCoef).^2)/sum(abs(s(k).recNScoef).^2);
end

Cg_dB = 10*log10(Cglobal);

% % Visualize
% numPointsPerCircle = size(s, 2); numCircles = size(s, 3);
% visualObj = animation({1:size(s, 1), 1:numPointsPerCircle, 1:numCircles},...
%     {Cg_dB}, {'Type', 'Points', 'Circle'}, {'Cancellation'}, [], []);

% Generate graph for TFM report
ax = axes(figure);
Cg_dB_scal = permute(Cg_dB(2, :, :), [2 3, 1]);
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
% widthInPixels = 600;
% heightInPixels = 600;
% fig = ax.Parent;
% fig.Position(3:4) = [widthInPixels, heightInPixels];
% fig.Position(1:2) = [0 0];

fontSize_axesWidth_ratio = 0.08;
fontSize = ax.Position(3) * fontSize_axesWidth_ratio;
ax.XLabel.FontUnits = 'normalized';
ax.XLabel.FontSize = fontSize;
ax.YLabel.FontUnits = 'normalized';
ax.YLabel.FontSize = fontSize;

printfig(ax.Parent, imagesPath, 'Experiment1_globalCancDifNSpos', 'eps');
%% Visualization: correction factor AllTogether
corrFact1 = zeros(numPointsPerCircle, numCircles);
for o = 1:numCircles
    for p = 1:numPointsPerCircle
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

[meanPhase, ~, stdPhase] = circularDistributionParameters(angle(data));
stdPhase = rad2deg(stdPhase);
meanPhase = rad2deg(meanPhase);
normStdPhase = stdPhase./meanPhase;

str = cell(size(data, 2), 1);
for k = 1:size(data, 2)
    str{k} = sprintf('$|\\overline{x}| = %.2g$, angle$(\\overline{x}) = %.0f$, $\\sigma_{n} = %.2g$, $\\sigma_{n, abs} = %.2g$, $\\sigma_{n, phase} = %.2g$', meanAbs(k), meanPhase(k), normStd(k), normStdAbs(k), normStdPhase(k));
end
l = legend(ax, str);
l.Interpreter = 'latex';

ax = axes(figure);
yyaxis(ax, 'left')
plot(ax, rad2deg([alpha; 2*pi]), abs([data; data(1,:)]))
ax.XLabel.String = '\alpha';
ax.YLim = [0, max(abs(data(:)))*2];
ax.XLim = [0, 360];
ax.YLabel.String = '|\Psi|';
yyaxis(ax, 'right');
lAngle = plot(ax, rad2deg([alpha; 2*pi]), rad2deg(unwrap(angle([data; data(1,:)]))));
ax.YLim = [-180 180];
ax.YLabel.String = 'Phase of \Psi (º)';

legLab = cell(numCircles, 1);
for c = 1:numCircles
    legLab{c} = num2str(radius(c));
end
l = legend(ax, [legLab; legLab]);
l.Interpreter = 'latex';
title(l, 'R (m)')

fontSize_axesWidth_ratio = 0.08;
fontSize = ax.Position(3) * fontSize_axesWidth_ratio;
ax.XLabel.FontUnits = 'normalized';
ax.XLabel.FontSize = fontSize;
ax.YAxis(1).Label.FontUnits = 'normalized';
ax.YAxis(1).Label.FontSize = fontSize;
ax.YAxis(2).Label.FontUnits = 'normalized';
ax.YAxis(2).Label.FontSize = fontSize;

printfig(ax.Parent, imagesPath, 'Experiment1_globalCancScaleFactor', 'eps');