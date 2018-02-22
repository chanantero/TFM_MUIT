function svgText = WFSarraySVG(varargin)

p = inputParser;

addRequired(p, 'viewBox');
addParameter(p, 'NSposition', [0 0]);
addParameter(p, 'NSangle', 0);
addParameter(p, 'microPositions', double.empty(0, 2));
addParameter(p, 'WFSdiagramIndex', []);
addParameter(p, 'WFSdiagramName', []);
addParameter(p, 'microDiagramIndex', []);
addParameter(p, 'microDiagramName', []);
addParameter(p, 'backgroundFileName', []);

parse(p, varargin{:});

viewBox = p.Results.viewBox;
NSpos = p.Results.NSposition;
NSangle = p.Results.NSangle;
microPos = p.Results.microPositions;
WFSdiagInd = p.Results.WFSdiagramIndex;
WFSdiagName = p.Results.WFSdiagramName;
microDiagInd = p.Results.microDiagramIndex;
microDiagName = p.Results.microDiagramName;
backgroundName = p.Results.backgroundFileName;

if ischar(WFSdiagName)
    WFSdiagName = {WFSdiagName};
end

if ischar(microDiagName)
    microDiagName = {microDiagName};
end

% Generate positions and orientations of WFS array loudspeakers
d = 0.18; % Separation between two contiguous loudspeakers. Size of one loudspeaker
nb = 8; % Bottom and upper sides of the octogon (2 sides)
nd = 8; % Diagonal sides of the octogon (4 sides)
nl = 24; % Lateral side of the octogon (2 sides)
betabd = 45; % Deviation angle between bottom/upper and diagonal sides
[ xLoud, yLoud, angleLoud ] = octogon(d, nb, nd, nl, betabd);

% Position and orientation of noise loudspeaker
xNoise = NSpos(:, 1); yNoise = NSpos(:, 2); angleNoise = NSangle; noiseName = 'loudspeakerSound';

% Positions of microphones
xMicro = microPos(:, 1);
yMicro = microPos(:, 2);

% Position and orientation of radiation diagrams
xDiag = [xLoud(WFSdiagInd); xMicro(microDiagInd)]; yDiag = [yLoud(WFSdiagInd), yMicro(microDiagInd)]; angleDiag = [angleLoud(6); 0];
diagName = [WFSdiagName; microDiagName];

% Viewbox and viewport
viewPort2viewBoxRatio = 100;
width = viewBox(3)*viewPort2viewBoxRatio;
height = viewBox(4)*viewPort2viewBoxRatio;

% Create string that will draw them in SVG
drawLoudspeakerStr = '<use xlink:href="symbolLibrary.svg#%s" x="%g" y="%g" transform="rotate(%g %g %g)"/>\n';
drawMicrophoneStr = '<use xlink:href="symbolLibrary.svg#microphone" transform="translate(%g %g) scale(%g)"/>\n';
drawDiagramStr = '<use xlink:href="symbolLibrary.svg#%s" x="%g" y="%g" transform="rotate(%g %g %g)" />\n';
drawBackgroundStr = '<image xlink:href="%s" x="%g" y="%g" width="%g" height="%g" preserveAspectRatio="none"/>';

if ~isempty(backgroundName)
    strBackground = sprintf(drawBackgroundStr, backgroundName, viewBox(1), viewBox(2), viewBox(3), viewBox(4));
else
    strBackground = '';
end

strLoud = '';
numLoudspeakers = numel(xLoud);
for k = 1:numLoudspeakers
    strLoud = [strLoud, sprintf(drawLoudspeakerStr, 'loudspeaker', xLoud(k), yLoud(k), angleLoud(k), xLoud(k), yLoud(k))];
end
strLoud = [sprintf('<g>\n'), strLoud, sprintf('</g>\n')];

strLoudNoise = '';
numNoiseLoudspeakers = numel(xNoise);
for k = 1:numNoiseLoudspeakers
    strLoudNoise = [strLoudNoise, sprintf(drawLoudspeakerStr, noiseName, xNoise(k), yNoise(k), angleNoise(k), xNoise(k), yNoise(k))];
end

strMicro = '';
numMicrophones = numel(xMicro);
for k = 1:numMicrophones
    strMicro = [strMicro, sprintf(drawMicrophoneStr, xMicro(k), yMicro(k), 0.8)];
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

% Insert background image
svgText = strrep(svgText, '[background]', strBackground);

% Insert loudspeakers in the desired positions and orientations
svgText = strrep(svgText, '[WFSarray]', strLoud);

% Noise Loudspeakers
svgText = strrep(svgText, '[NoiseLouds]', strLoudNoise);

% Insert string to draw microphones
svgText = strrep(svgText, '[MicrophoneArray]', strMicro);

% Insert diagrams
svgText = strrep(svgText, '[RadDiag]', strDiag);

end

