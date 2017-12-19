%% Draw WFS array in SVG

% Generate positions and orientations of WFS array loudspeakers
d = 0.18; % Separation between two contiguous loudspeakers. Size of one loudspeaker
nb = 8; % Bottom and upper sides of the octogon (2 sides)
nd = 8; % Diagonal sides of the octogon (4 sides)
nl = 24; % Lateral side of the octogon (2 sides)
betabd = 45; % Deviation angle between bottom/upper and diagonal sides
[ xLoud, yLoud, angleLoud ] = octogon(d, nb, nd, nl, betabd);

% Position and orientation of noise loudspeaker
xNoise = -1; yNoise = 2; angleNoise = 0; noiseName = 'loudspeakerSound';

% Positions of microphones
xMicro = [2; 2];
yMicro = [1; 5];

% Position and orientation of radiation diagrams
xDiag = [xLoud(6); xMicro(2)]; yDiag = [yLoud(6), yMicro(2)]; angleDiag = [angleLoud(6); 0];
diagName = {'diagramDipole', 'diagramMonopoleNoise'};

% Viewbox and viewport
width = 600; height = 800;
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
for k = 1:numLoudspeakers
    strLoud = [strLoud, sprintf(drawLoudspeakerStr, 'loudspeaker', xLoud(k), yLoud(k), angleLoud(k), xLoud(k), yLoud(k))];
end
strLoud = [sprintf('<g>\n'), strLoud, sprintf('</g>\n')];

strLoudNoise =  sprintf(drawLoudspeakerStr, noiseName, xNoise, yNoise, angleNoise, xNoise, yNoise);

strMicro = '';
numMicrophones = numel(xMicro);
for k = 1:numMicrophones
    strMicro = [strMicro, sprintf(drawMicrophoneStr, xMicro(k), yMicro(k))];
end

strDiag = '';
numDiagrams = numel(xDiag);
for k = 1:numDiagrams
    strDiag = [strDiag, sprintf(drawDiagramStr, diagName{k}, xDiag(k), yDiag(k), angleDiag(k), xDiag(k), yDiag(k))];
end

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

% Insert string to draw microphones
svgText = strrep(svgText, '[MicrophoneArray]', strMicro);

% Insert diagrams
svgText = strrep(svgText, '[RadDiag]', strDiag);

% Write file
path = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\';
name = 'WFSarrayScheme';
destFile = fopen([path, name, '.svg'], 'w', 'n', 'UTF-8');
fwrite(destFile, svgText, 'char');
fclose(destFile);

% Export to PNG and EMF
% Matlab should be in the folder with the .svg to avoid path problems
system(['inkscape -z "', path, fileName, '.svg"', ' -e "', path, fileName, '.png"', ' --export-dpi=384']);
system(['inkscape -z "', path, fileName, '.svg"', ' --export-emf "', path, fileName, '.emf"']);