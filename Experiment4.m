%% Experiment 4
% Noise source outside the plane of the WFS array. Change z coordinate

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

%% Apply cancellation

% Circles
numPointsPerCircle = 1;
radius = [5 50 5000];
zAngleLim = [0 pi/4]; % Angle between the WFS plan and the vector from the center of the array to the noise source
numPointsZ = 25;
numCircles = numel(radius);
alpha = linspace(0, 2*pi, numPointsPerCircle + 1); alpha = alpha(1:end-1)';
zAngleVec = linspace(zAngleLim(1), zAngleLim(2), numPointsZ);
[Alpha, Radius, Zangle] = ndgrid(alpha, radius, zAngleVec);
x = centreX + Radius(:).*cos(Alpha(:));
y = centreY + Radius(:).*sin(Alpha(:));
z = Radius(:).*tan(Zangle(:));

NSpositions = [x, y, z];
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
s = reshape(s, [2, numPointsPerCircle, numCircles, numPointsZ]);

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
% visualObj = animation({1:size(s, 1), 1:numPointsPerCircle, 1:numCircles, 1:numPointsZ},...
%     {Cg_dB}, {'Type', 'Points', 'Circle', 'Z'}, {'Cancellation'}, [], []);

% Generate graph for TFM report
ax = axes(figure);
Cg_dB_scal = permute(Cg_dB(2, 1, :, :), [4, 3, 2, 1]);
plot(ax, rad2deg(zAngleVec), Cg_dB_scal)
ax.XLabel.String = '\theta (º)';
ax.YLabel.String = 'Global cancellation (dB)';
legLab = cell(numCircles, 1);
for c = 1:numCircles
    legLab{c} = num2str(radius(c));
end
l = legend(ax, legLab);
title(l, 'R')

% % Print graph
% fontSize_axesWidth_ratio = 0.08;
% fontSize = ax.Position(3) * fontSize_axesWidth_ratio;
% ax.XLabel.FontUnits = 'normalized';
% ax.XLabel.FontSize = fontSize;
% ax.YLabel.FontUnits = 'normalized';
% ax.YLabel.FontSize = fontSize;
% 
% printfig(ax.Parent, imagesPath, 'Experiment4_globalCancDifNSz', 'eps');
%% Visualization: correction factor AllTogether
sizeS = size(s);
corrFact = zeros(sizeS(2:end));
for k = 1:prod(sizeS(2:end))
    % Find WFS coefficient that is not 0
    defaultWFScoef = s(1, k).WFScoef;
    ind = find(defaultWFScoef ~= 0, 1, 'first');
    % Calculate
    corrFact(k) = mean(s(2, k).WFScoef(ind)./defaultWFScoef(ind));
end

ax = axes(figure);
ax.Title.Interpreter = 'latex';
hold on
data = corrFact;
cmap = colormap('lines');
data = mergeAndPermute(data, {3, [2 1]});
xData = real(data);
yData = imag(data);
for k = 1:numCircles
    xData = real(data(:, k));
    yData = imag(data(:, k));
    scat = scatter(ax, xData, yData, 10, cmap(k,:), 'filled');
end
maxAbs = max(abs(data(:)));
ax.XLim = [-maxAbs, maxAbs]; ax.YLim = [-maxAbs, maxAbs];
ax.DataAspectRatio = [1 1 3];

%% Image of different NS angles with matlab

f = figure;
ax = copyobj(obj.ax, f);
ax.CameraPosition =  [21.4837  -39.5552   21.3544];
ax.NextPlot = 'Add';

scatRec = findobj(ax, 'Tag', 'receiver');
delete(scatRec);

scatWFS = findobj(ax, 'Tag', 'loudspeakers');
scatWFS.CData = [0 0 0];

NSposition = [5 7 1];
scatNS = findobj(ax, 'Tag', 'source');
scatNS.MarkerFaceColor = 'flat';
scatNS.XData = NSposition(1);
scatNS.YData = NSposition(2);
scatNS.ZData = NSposition(3);
scatNS.CData = [0 0 0];

plot3(ax, [centerX, NSposition(1)], [centerY, NSposition(2)], [0, NSposition(3)], '--k')
plot3(ax, [centerX, NSposition(1)], [centerY, NSposition(2)], [0, 0], '--k')

maxX = max([NSposition(1); scatWFS.XData(:)]);
minX = min([NSposition(1); scatWFS.XData(:)]);
ax.XLim = [minX - 0.1, maxX + 0.1];

maxY = max([NSposition(2); scatWFS.YData(:)]);
minY = min([NSposition(2); scatWFS.YData(:)]);
ax.YLim = [minY - 0.1, maxY + 0.1];

maxZ = max([NSposition(3); scatWFS.ZData(:)]);
minZ = min([NSposition(3); scatWFS.ZData(:)]);
ax.ZLim = [minZ - 0.1, maxZ + 0.1];

ax.Box = 'off';
ax.XColor = 'none';
ax.YColor = 'none';
ax.ZColor = 'none';

%% Image of different NS angles with SVG

WFSpositions = obj.WFSposition;

extRectXmin = min(obj.WFSposition(:, 1));
extRectXmax = max(obj.WFSposition(:, 1));
extRectYmin = min(obj.WFSposition(:, 2));
extRectYmax = max(obj.WFSposition(:, 2));
octagonRectPos = [extRectXmin, extRectYmin, extRectXmax - extRectXmin, extRectYmax - extRectYmin];
centerX = (extRectXmax + extRectXmin)/2;
centerY = (extRectYmax + extRectYmin)/2;
centre = [centerX, centerY];

theta = 60;
phi = -60;
alpha = phi + 90;

NSposition = [5 7 1];

% Consult appendix of thesis
h = 100;
ro = h*[cosd(phi)*sind(theta), sind(phi)*sind(theta), cosd(theta)] + [centerX, centerY, 0];
l = [sind(phi); -cosd(phi); 0];
u = [-cosd(phi)*cosd(theta); -sind(phi)*cosd(theta); sind(theta)];

rs_NS = [-l, u]'*(NSposition)' + centre';

% rs_centre = [-l, u]'*centre';
% strLine1 = makePath('000000', 0.01, centreX, centreY, NSpositions(selLoud,1), NSpositions(selLoud,2), 'line');


margin = 3;
viewBox = [octagonRectPos(1) - margin, octagonRectPos(2) - margin, octagonRectPos(3) + margin*2, octagonRectPos(4) + margin*2];
objSVG = SVGdrawer('viewBox', viewBox, 'NSpositions', rs_NS', 'microSymbol', 'dot', 'microSize', 0.01);

svgText = objSVG.getSVG;


% Mannually apply transformation to loudspeaker group

translate = (centerY - extRectYmin)*(1 - cosd(theta));
strTransf_temp = 'transform="translate(0 %g) scale(1 %g) rotate(%g %g %g)"';
strTransf = sprintf(strTransf_temp, translate, cosd(theta), -alpha, centerX, centerY);

tag = 'id="WFSarray"';
k = strfind(svgText, tag);
svgText = [svgText(1:k+numel(tag)-1), ' ', strTransf, svgText(k+numel(tag):end)];

% Write file
name = 'Experiment4_diffNSzpos';
fileName = [imagesPath, name, '.svg'];
destFile = fopen(fileName, 'w', 'n', 'UTF-8');
fwrite(destFile, svgText, 'char');
fclose(destFile);

% Convert it to pdf
currentFolder = pwd;
cd(imagesPath); % Needed for inkscape to link svg files properly
system(['inkscape -z "', imagesPath, name, '.svg" --export-pdf="', imagesPath, name, '.pdf"'])
cd(currentFolder)