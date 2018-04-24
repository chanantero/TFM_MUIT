%% Experiment 5
% Impulse responses for theoretical case

%% Preamble
pathSetUp;

imagesPath = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\TFM\Img\';

dataPathName = [globalPath, 'Data\'];
ID = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

%% System set up.
obj = SimulationController;

obj.NSposition = [3.35 -0.2 0]; % Assumed real position
obj.amplitude = 1;
obj.amplitude(2) = -obj.amplitude(1);
obj.phase = 0;
obj.frequency = 440;
obj.Fs = 44100;

load('WFSTool/WFSfilter.mat')
obj.WFSToolObj.freqFilter = hTotal;

% Microphone positions
% Rectangular grid
marginRatio = 0.6;
numPointsX = 10;
numPoinstY = 10;
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
obj.WFSToolObj.domain = 'time';
obj.setAcousticPaths('NS', 'theoretical', 'WFS', 'theoretical');

%% Signal to be transmitted by the noise source
f = 440;
t = 0:1/obj.WFSToolObj.Fs:2;
x = cos(2*pi*f*t);

obj.NScoef = x;
obj.NSVcoef = -x;

obj.WFSToolObj.WFScalculation();
