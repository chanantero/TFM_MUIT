%% Draw WFS calculation parameters in SVG

% Generate positions and orientations of WFS array loudspeakers
d = 0.18; % Separation between two contiguous loudspeakers. Size of one loudspeaker
nb = 8; % Bottom and upper sides of the octogon (2 sides)
nd = 8; % Diagonal sides of the octogon (4 sides)
nl = 24; % Lateral side of the octogon (2 sides)
betabd = 45; % Deviation angle between bottom/upper and diagonal sides
[ xLoud, yLoud, angleLoud ] = octogon(d, nb, nd, nl, betabd);

% Select only the loudspeakers from the low left corner
ind = [1:16, 89:96];
xLoud = xLoud(ind);
yLoud = yLoud(ind);
angleLoud = angleLoud(ind);

% Position and orientation of noise loudspeaker
xNoise = -0.6; yNoise = 0.7; angleNoise = 0; noiseName = 'loudspeakerSound';

% Position and orientation of radiation diagrams
selLoud = 3; % Index of selected loudspeaker
xDiag = xLoud(selLoud); yDiag = yLoud(selLoud); angleDiag = angleLoud(selLoud);
diagName = {'diagramDipole'};

% Viewbox and viewport
width = 400; height = 400;
viewPort2viewBoxRatio = 100;

viewBoxWidth = width/viewPort2viewBoxRatio;
viewBoxHeight = height/viewPort2viewBoxRatio;

rangeX = max(xLoud) - min(xLoud);
rangeY = max(yLoud) - min(yLoud);

marginX = (viewBoxWidth - rangeX)/2;
marginY = (viewBoxHeight - rangeY)/2;

viewBox = [min(xLoud) - marginX, min(yLoud) - marginY, viewBoxWidth, viewBoxHeight];

% Create string that will draw them in SVG
drawLoudspeakerStr = '<use xlink:href="symbolLibrary.svg#%s" x="%g" y="%g" transform="rotate(%g %g %g)"/>\n';
drawMicrophoneStr = '<use xlink:href="symbolLibrary.svg#microphone" x="%g" y="%g"/>\n';
drawDiagramStr = '<use xlink:href="symbolLibrary.svg#%s" x="%g" y="%g" transform="rotate(%g %g %g)" />\n';

strLoud = '';
numLoudspeakers = numel(xLoud);
for k = 1:16
    strLoud = [strLoud, sprintf(drawLoudspeakerStr, 'loudspeakerSound', xLoud(k), yLoud(k), angleLoud(k), xLoud(k), yLoud(k))];
end
for k = 17:numLoudspeakers
    strLoud = [strLoud, sprintf(drawLoudspeakerStr, 'loudspeaker', xLoud(k), yLoud(k), angleLoud(k), xLoud(k), yLoud(k))];
end
strLoud = [sprintf('<g>\n'), strLoud, sprintf('</g>\n')];

strLoudNoise =  sprintf(drawLoudspeakerStr, noiseName, xNoise, yNoise, angleNoise, xNoise, yNoise);

strDiag = '';
numDiagrams = numel(xDiag);
for k = 1:numDiagrams
    strDiag = [strDiag, sprintf(drawDiagramStr, diagName{k}, xDiag(k), yDiag(k), angleDiag(k), xDiag(k), yDiag(k))];
end

% Draw line from noise loudspeaker to selected WFS array loudspeaker. Make
% it longer to mark the alpha angle
% d = norm([xLoud(selLoud) - xNoise, yLoud(selLoud) - yNoise], 2); % Distance between loudspeakers
l = 1; % Additional length of line
dir = [xLoud(selLoud) - xNoise, yLoud(selLoud) - yNoise];
d = norm(dir); % Distance
endPoint = [xNoise, yNoise] + dir + dir/d*l;
strLine = makePath('000000', 0.01, xNoise, yNoise, endPoint(1), endPoint(2), 'line');

% Draw curly brace to indicate the distance
origWidth = 52.179153; % Original width of the curly brace
origHeight = 163.10518; % Original height of the curly brace
xEnd = xNoise; yEnd = yNoise; xStart = xLoud(selLoud); yStart = yLoud(selLoud);
v = [xEnd - xStart; yEnd - yStart];
L = norm(v); % Length
scaleFactor = L/origHeight;
alpha = atan2d(v(2), v(1));
angleRot = alpha - 90;
margin = 0.15/scaleFactor; % Separation of curly brace from the thing is measuring
strCurlyBrace = sprintf('<use xlink:href="symbolLibrary.svg#leftBracket" transform="translate(%g %g) rotate(%g 0 0) scale(%g) translate(%g 0)"/>',...
    xStart, yStart, angleRot, scaleFactor, -(origWidth + margin));

% Draw the symbol of the distance
origPeak = [-0.15/scaleFactor origHeight/2]; % Position of peak of the curly brace
SVGcommand = sprintf('translate(%g %g) rotate(%g 0 0) scale(%g) translate(%g 0)',...
    xStart, yStart, angleRot, scaleFactor, -(origWidth + margin));
transfMatrix = transformationMatrix(SVGcommand);
newPoint = transfMatrix * [origPeak'; 1];
strDistance = makeText( '$d$', 0.1, '000000', newPoint(1), newPoint(2), 'distanceSymbol');

% Draw line to indicate broadside direction
l = 1; % Length of line
xEnd = xLoud(selLoud) + cosd(angleLoud(selLoud))*l;
yEnd = yLoud(selLoud) + sind(angleLoud(selLoud))*l;
strBroadside = makePath('000000', 0.01, xLoud(selLoud), yLoud(selLoud), xEnd, yEnd, 'broadside');

% Draw arc to mark the angle alpha
endAngle = angleLoud(selLoud);
startAngle = rad2deg(cart2pol(dir(1), dir(2)));
strArc = makeCircumferenceArc(xLoud(selLoud), yLoud(selLoud), 0.25, startAngle, endAngle, '000000', 0.01);

% Write the symbols of the angle and the distance
strAlpha = makeText( '$\alpha$', 0.1, '000000', xLoud(selLoud) + 0.3, yLoud(selLoud) + 0.15, 'alphaSymbol');

% Read file
file = fopen('C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\WFSarrayScheme_Template.svg');
svgText = fread(file, Inf, '*char')';
fclose(file);

% Set viewport and viewbox
svgText = strrep(svgText, '[viewPortWidth]', num2str(width));
svgText = strrep(svgText, '[viewPortHeight]', num2str(height));
svgText = strrep(svgText, '[viewBox]', num2str(viewBox));

% Insert loudspeakers in the desired positions and orientations
svgText = strrep(svgText, '[WFSarray]', strLoud);

% Noise Loudspeakers
svgText = strrep(svgText, '[NoiseLouds]', strLoudNoise);

% Insert diagrams
svgText = strrep(svgText, '[RadDiag]', strDiag);

% Delete prewritten block for microphones
svgText = strrep(svgText, '[MicrophoneArray]', '');

% Insert lines, arc and text
svgText = strrep(svgText, '[Other]', [strLine, sprintf('\n'), strBroadside, ...
    sprintf('\n'), strArc, sprintf('\n'), strAlpha, sprintf('\n'), strCurlyBrace, ...
    sprintf('\n'), strDistance]);

% Write file
path = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\';
name = 'WFSparameters';
destFile = fopen([path, name, '.svg'], 'w', 'n', 'UTF-8');
fwrite(destFile, svgText, 'char');
fclose(destFile);

% % Export to PNG and EMF
% % Matlab should be in the folder with the .svg to avoid path problems
% system(['inkscape -z "', path, name, '.svg"', ' -e "', path, name, '.png"', ' --export-dpi=384']);
% system(['inkscape -z "', path, name, '.svg"', ' --export-emf "', path, name, '.emf"']);
system(['inkscape -z "', path, name, '.svg" --export-pdf="', path, name, '.pdf" --export-latex'])