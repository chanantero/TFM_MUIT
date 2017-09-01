%% Create axes
% Import image
[imWFSArray, map, transparency] = imread('C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\WFSarrayScheme.png', 'png');
[numRows, numColumns, numColors] = size(imWFSArray);

xLim = [viewBox(1), viewBox(1) + viewBox(3)];
yLim = [viewBox(2), viewBox(2) + viewBox(4)];

ax = axes(figure);
image(ax, 'XData', xLim, 'YData', yLim, 'CData', imWFSArray, 'AlphaData', transparency);
ax.DataAspectRatio = [1, 1, 1];
ax.XLim = xLim;
ax.YLim = yLim;
ax.XTick = [];
ax.YTick = [];

ax.CLim = [-0.5, 0.5];
basicColormap = [0 1 0; 1 1 1; 1 1 0];
colorMap_R = interp1(1:3, basicColormap(:,1), linspace(1, 3, 128)');
colorMap_G = interp1(1:3, basicColormap(:,2), linspace(1, 3, 128)');
colorMap_B = interp1(1:3, basicColormap(:,3), linspace(1, 3, 128)');
extendedColormap = [colorMap_R, colorMap_G, colorMap_B];
colormap(ax, extendedColormap);

%% Simulate first scenario: only real source
numLoudspeakers = 96;
positionNoise = [-1 2 0];
amplitude = 0.5;
phase = 0;
frequency = 600;

% Set scenario
obj = WFSToolSimple();
obj.changeScenario(96);

obj.setNumNoiseSources(1);
obj.noiseSourceChannelMapping = 1;

obj.amplitude = amplitude;
obj.phase = phase;
obj.frequency = frequency;
obj.noiseSourcePosition = positionNoise;
obj.setVirtual(false);
obj.setReal(true);
obj.updateReprodPanelBasedOnVariables();

obj.WFScalculation();

% Configure simulator
delete(obj.simulObj.imag);
obj.simulObj.ax = ax;
obj.simulObj.XLim = xLim;
obj.simulObj.YLim = yLim;
obj.simulObj.XnumPoints = round(numColumns/8);
obj.simulObj.YnumPoints = round(numRows/8);
obj.simulObj.generateMeasurePoints();

obj.simulObj.simulate();

% Save BMP
indices = scaled2indexedColors(size(extendedColormap, 1), ax.CLim, obj.simulObj.imag.CData);
imwrite(indices, extendedColormap, 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\wave.bmp');

% Generate GIF
numFramesPerPeriod = 10; freq = 600;
t = (0:numFramesPerPeriod-1)/(freq*numFramesPerPeriod);
numFrames = numel(t);
imagGIF = zeros(obj.simulObj.YnumPoints, obj.simulObj.XnumPoints, 1, numFrames);
for k = 1:numFrames
U = obj.simulObj.field.*repmat(exp(1i*2*pi*obj.simulObj.freq'*t(k)), [obj.simulObj.numMeasurePoints, 1]);
U = reshape(sum(U, 2), obj.simulObj.XnumPoints, obj.simulObj.YnumPoints).';
imagGIF(:, :, 1, k) = real(U);
end
indices = scaled2indexedColors(size(extendedColormap, 1), ax.CLim, imagGIF);

secondsPerPeriod = 1;

imwrite(indices, extendedColormap, 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\wave.gif', ...
    'LoopCount', 100, 'DelayTime', secondsPerPeriod/numFramesPerPeriod);

%% Simulate second scenario: only virtual noise and alltogether
obj.setVirtual(true);
obj.setReal(true);
obj.updateReprodPanelBasedOnVariables();
obj.WFScalculation();
obj.WFS2realRatio();

obj.simulObj.simulate();

path = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\';
fileName = 'waveCancellation';
% Save BMP
indices = scaled2indexedColors(size(extendedColormap, 1), ax.CLim, obj.simulObj.imag.CData);
imwrite(indices, extendedColormap, [path, fileName, '.bmp']);

% Generate GIF
numFramesPerPeriod = 30; freq = 600;
t = (0:numFramesPerPeriod-1)/(freq*numFramesPerPeriod);
numFrames = numel(t);
imagGIF = zeros(obj.simulObj.YnumPoints, obj.simulObj.XnumPoints, 1, numFrames);
for k = 1:numFrames
U = obj.simulObj.field.*repmat(exp(1i*2*pi*obj.simulObj.freq'*t(k)), [obj.simulObj.numMeasurePoints, 1]);
U = reshape(sum(U, 2), obj.simulObj.XnumPoints, obj.simulObj.YnumPoints).';
imagGIF(:, :, 1, k) = real(U);
end
indices = scaled2indexedColors(size(extendedColormap, 1), ax.CLim, imagGIF);

secondsPerPeriod = 1;

imwrite(indices, extendedColormap, [path, fileName, '.gif'], ...
    'LoopCount', 100, 'DelayTime', secondsPerPeriod/numFramesPerPeriod);

% Only virtual
realChannels = obj.noiseSourceChannelMapping(obj.real);
realFlag = ismember(obj.loudspeakerChannelMapping, realChannels);
WFSflag = ~realFlag;
obj.loudspeakerCoefficient(realFlag, :) = 0;

obj.simulObj.simulate();

path = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\';
fileName = 'waveWFS';

% Save BMP
indices = scaled2indexedColors(size(extendedColormap, 1), ax.CLim, obj.simulObj.imag.CData);
imwrite(indices, extendedColormap, [path, fileName, '.bmp']);

% Generate GIF
numFramesPerPeriod = 30; freq = 600;
t = (0:numFramesPerPeriod-1)/(freq*numFramesPerPeriod);
numFrames = numel(t);
imagGIF = zeros(obj.simulObj.YnumPoints, obj.simulObj.XnumPoints, 1, numFrames);
for k = 1:numFrames
U = obj.simulObj.field.*repmat(exp(1i*2*pi*obj.simulObj.freq'*t(k)), [obj.simulObj.numMeasurePoints, 1]);
U = reshape(sum(U, 2), obj.simulObj.XnumPoints, obj.simulObj.YnumPoints).';
imagGIF(:, :, 1, k) = real(U);
end
indices = scaled2indexedColors(size(extendedColormap, 1), ax.CLim, imagGIF);

secondsPerPeriod = 1;

imwrite(indices, extendedColormap, [path, fileName, '.gif'], ...
    'LoopCount', 100, 'DelayTime', secondsPerPeriod/numFramesPerPeriod);