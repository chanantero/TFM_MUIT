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
obj.frequency = 800;
obj.Fs = 44100;

load('WFSTool/WFSfilter.mat')
obj.WFSToolObj.freqFilter = hTotal;

% Microphone positions
% Rectangular grid
marginRatio = 0.6;
numPointsX = 2;
numPoinstY = 2;
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


%% Signal to be transmitted by the noise source
f = obj.frequency;
t = 0:1/obj.WFSToolObj.Fs:1;
x = obj.amplitude(1) * cos(2*pi*f*t + obj.phase(1));
x(t > 1) = 0;

obj.WFSToolObj.domain = 'time';
obj.NScoef = x;
obj.NSVcoef = -x;

%% Simulate
% Acoustic paths
obj.setAcousticPaths('NS', 'theoretical', 'WFS', 'theoretical');

obj.WFSToolObj.updateFiltersWFS();

% Simulate only the noise source
obj.WFSToolObj.virtual = [false; false];
obj.WFSToolObj.WFScalculation();
obj.WFSToolObj.simulate();
simulFieldOnlyNoise = obj.WFSToolObj.simulField;

% Simulate all together
obj.WFSToolObj.virtual = [false; true];
obj.WFSToolObj.WFScalculation();
obj.WFSToolObj.simulate();
simulField = obj.WFSToolObj.simulField;

% save([dataPathName, 'pruebaTiempo', ID, '.mat'], 'simulFieldOnlyNoise', 'simulField')

% % Visualize
% ax = axes(figure, 'NextPlot', 'Add');
% lOnlyNoise = plot(ax, simulFieldOnlyNoise(3, :));
% lAll = plot(ax, simulField(1, :));
% lDiff = plot(ax, simulField(1, :) - simulFieldOnlyNoise(1, :));
% 
% ax.Children = flip(ax.Children);
% 
% indRec = 4;
% lOnlyNoise.YData = simulFieldOnlyNoise(indRec, :);
% lAll.YData = simulField(indRec, :);
% lDiff.YData = simulField(indRec, :) - simulFieldOnlyNoise(indRec, :);
% 
% % plot(ax, simulField(1,:)./simulFieldOnlyNoise(1,:))
% plot(ax, simulField(1,:)./simulFieldOnlyNoise(1,:))

%% Identify IQ component
recNScoef_time = signal2pulseCoefficientMatrix([0 1], f, 1, simulFieldOnlyNoise', obj.WFSToolObj.Fs);
recWFScoef_time = signal2pulseCoefficientMatrix([0.1 0.9], f, 1, (simulField - simulFieldOnlyNoise)', obj.WFSToolObj.Fs);
recCoef_time = signal2pulseCoefficientMatrix([0.1 0.9], f, 1, simulField', obj.WFSToolObj.Fs);
WFScoef_time = signal2pulseCoefficientMatrix([0.1 0.9], f, 1, obj.WFSToolObj.WFSarrayCoefficient(81:96, :)', obj.WFSToolObj.Fs);
NScoef_time = signal2pulseCoefficientMatrix([0.1 0.9], f, 1, x', obj.WFSToolObj.Fs);

s_time = SimulationController.generateExportStructure(...
                'NSRcoef', NScoef_time,...
                'NSVcoef', -NScoef_time,...
                'WFScoef', WFScoef_time,...
                'microCoef', recCoef_time,...
                'microCoefNS', recNScoef_time,...
                'microCoefWFS', recWFScoef_time,...
                'NSRpos', obj.NSRposition,...
                'NSVpos', obj.NSVposition,...
                'WFSpos', obj.WFSposition,...
                'microPos', obj.microPos,...
                'Frequency', obj.frequency...
                );

%% Result with frequency responses
obj.WFSToolObj.domain = 'frequency';
obj.setAcousticPaths('NS', 'theoretical', 'WFS', 'theoretical');

obj.cancelResults = [];

% Simulate only the noise source
obj.WFSToolObj.virtual = [false; false];
obj.WFSToolObj.WFScalculation();
obj.cancel();

% Simulate all together
obj.WFSToolObj.virtual = [false; true];
obj.WFSToolObj.WFScalculation();
obj.cancel();

s_frequency = obj.cancelResults;

%% Compare frequency and time responses
abs(WFScoef_time')./abs(s(2).WFScoef(81:96))
abs(recNScoef_time')./abs(s(2).recNScoef)
abs(recCoef_time')./abs(s(2).recCoef)
abs(recWFScoef_time')./abs(s(2).recWFScoef)

angle(WFScoef_time.'./s(2).WFScoef(81:96))
angle(recNScoef_time.'./s(2).recNScoef)
angle(recCoef_time.'./s(2).recCoef)
angle(recWFScoef_time.'./s(2).recWFScoef)

abs((recWFScoef_time + recNScoef_time).')./abs(s(2).recCoef)
angle((recWFScoef_time + recNScoef_time).'./s(2).recCoef)

abs((recWFScoef_time + recNScoef_time).'./recCoef_time.')
angle((recWFScoef_time + recNScoef_time).'./recCoef_time.')

%% Visualization: 2D map, case by case

s = [s_frequency; s_time]; 

% Format structure
for p = 1:numel(s)
    s(p).NScoef = [s(p).NSRcoef; s(p).NSVcoef];
    s(p).NSposition = [s(p).NSRposition; s(p).NSVposition];
end

% Create simulationViewer object
objVis = simulationViewer(obj.ax, s);

%% Try it with a song
filename = 'C:\Users\Rubén\Music\Salsa\Gente de Zona - Traidora (Salsa Version)[Cover Audio] ft. Marc Anthony.mp3';
[y, Fs] = audioread(filename);
numSamp = size(y, 1);

% Filter to avoid aliasing
Y = fft(y(:, 1));
df = 1/numSamp;
freqs = 0:Fs/numSamp:(1-1/numSamp)*Fs;

fcut = 1800;
Y(freqs > fcut & freqs < Fs - fcut) = 0;

yFilt = ifft(Y);
% sound(yFilt, Fs)

% Crop it
tStart = 94; % seconds
dur = 6; % seconds
indStart = floor(tStart*Fs);
indEnd = ceil((tStart + dur)*Fs) - 1;
x = yFilt(indStart:indEnd)';
% sound(x, Fs)

% Assign signal
obj.NScoef = x;
obj.NSVcoef = -x;

% Simulate in steps so we don't run out of memory
    % Parameters
signalLength = length(x);
frameLength = 44100;
numFrames = ceil(signalLength/frameLength);
fillSize = numFrames*frameLength - length(x);
x = [x, zeros(1, fillSize)];

filterLength = length(obj.WFSToolObj.filtersWFS_IR);
[~, filterDelay] = max(obj.WFSToolObj.freqFilter);
filterDelay = filterDelay - 1;

acPathFilterLength = size(obj.WFSToolObj.noiseSourceAcousticPath, 3);

    % Initialize
obj.WFSToolObj.updateFiltersWFS();

previousOnlyNoiseWFSarrayCoef = zeros(obj.numWFS, filterLength - 1);
previousOnlyNoiseNoiseSourceCoef = zeros(obj.numNS, filterLength - 1);
previousWFSarrayCoef = zeros(obj.numWFS, filterLength - 1);
previousNoiseSourceCoef = zeros(obj.numNS, filterLength - 1);

previousOnlyNoise = zeros(obj.numMicro, acPathFilterLength - 1);
previous = zeros(obj.numMicro, acPathFilterLength - 1);

fieldOnlyNoise = zeros(obj.numMicro, frameLength*numFrames);
field = zeros(obj.numMicro, frameLength*numFrames);
    
    % Loop
for fr = 1:numFrames
    fprintf('%d/%d\n', fr, numFrames);
    
    % Compensate for the delay compensation performed in WFSToolSimple
    x_frame = [zeros(1, filterDelay),  x((fr-1)*frameLength + 1:frameLength*fr)];
    
    obj.NScoef = x_frame;
    obj.NSVcoef = -x_frame;

        % Simulate only the noise source
    obj.WFSToolObj.virtual = [false; false];
    obj.WFSToolObj.WFScalculation();
    
    WFSarraySignals = obj.WFSToolObj.WFSarrayCoefficient;
    noiseSourceSignals = obj.WFSToolObj.noiseSourceCoefficient_complete;
    
    WFSarraySignals(:, 1:filterLength - 1) = WFSarraySignals(:, 1:filterLength - 1) + previousOnlyNoiseWFSarrayCoef;
    previousOnlyNoiseWFSarrayCoef = WFSarraySignals(:, frameLength + 1:end);
    WFSarraySignals = WFSarraySignals(:, 1:frameLength);
    
    noiseSourceSignals(:, 1:filterLength - 1) = noiseSourceSignals(:, 1:filterLength - 1) + previousOnlyNoiseNoiseSourceCoef;
    previousOnlyNoiseNoiseSourceCoef = noiseSourceSignals(:, frameLength + 1:end);
    noiseSourceSignals = noiseSourceSignals(:, 1:frameLength);
    
    % Extend it
    WFSarraySignals = [WFSarraySignals, zeros(obj.numWFS, acPathFilterLength - 1)];
    noiseSourceSignals = [noiseSourceSignals, zeros(obj.numNS, acPathFilterLength - 1)];
    obj.WFSToolObj.WFSarrayCoefficient = WFSarraySignals;
    obj.WFSToolObj.noiseSourceCoefficient_complete = noiseSourceSignals;
    
    obj.WFSToolObj.simulate();
    fieldOnlyNoise_frame = obj.WFSToolObj.simulField;
    fieldOnlyNoise_frame(:, 1:acPathFilterLength - 1) = fieldOnlyNoise_frame(:, 1:acPathFilterLength - 1) + previousOnlyNoise;
    previousOnlyNoise = fieldOnlyNoise_frame(:, frameLength + 1:end);
    fieldOnlyNoise_frame = fieldOnlyNoise_frame(:, 1:frameLength);
    
        % Simulate all together
    obj.NScoef = x_frame;
    obj.NSVcoef = -x_frame;
    obj.WFSToolObj.virtual = [false; true];
    obj.WFSToolObj.WFScalculation();
    
    WFSarraySignals = obj.WFSToolObj.WFSarrayCoefficient;
    noiseSourceSignals = obj.WFSToolObj.noiseSourceCoefficient_complete;
    
    WFSarraySignals(:, 1:filterLength - 1) = WFSarraySignals(:, 1:filterLength - 1) + previousWFSarrayCoef;
    previousWFSarrayCoef = WFSarraySignals(:, frameLength + 1:end);
    WFSarraySignals = WFSarraySignals(:, 1:frameLength);
    
    noiseSourceSignals(:, 1:filterLength - 1) = noiseSourceSignals(:, 1:filterLength - 1) + previousNoiseSourceCoef;
    previousNoiseSourceCoef = noiseSourceSignals(:, frameLength + 1:end);
    noiseSourceSignals = noiseSourceSignals(:, 1:frameLength);
    
    % Extend it
    WFSarraySignals = [WFSarraySignals, zeros(obj.numWFS, acPathFilterLength - 1)];
    noiseSourceSignals = [noiseSourceSignals, zeros(obj.numNS, acPathFilterLength - 1)];
    obj.WFSToolObj.WFSarrayCoefficient = WFSarraySignals;
    obj.WFSToolObj.noiseSourceCoefficient_complete = noiseSourceSignals;
    
    obj.WFSToolObj.simulate();
    field_frame = obj.WFSToolObj.simulField;
    field_frame(:, 1:acPathFilterLength - 1) = field_frame(:, 1:acPathFilterLength - 1) + previous;
    previous = field_frame(:, frameLength + 1:end);
    field_frame = field_frame(:, 1:frameLength);
    
    fieldOnlyNoise(:, (fr-1)*frameLength + 1:frameLength*fr) = fieldOnlyNoise_frame;
    field(:, (fr-1)*frameLength + 1:frameLength*fr) = field_frame;

end

ax = axes(figure, 'NextPlot', 'Add');
% plot(ax, fieldOnlyNoise_frame(1, :))
% plot(ax, field_frame(1, :))
% plot(ax, noiseSourceSignals(2,:))
% plot(ax, WFSarraySignals(90, :))
lOnlyNoise = plot(ax, fieldOnlyNoise(1,:));
lAll = plot(ax, field(1,:));

indRec = 4;
lOnlyNoise.YData = fieldOnlyNoise(indRec, :);
lAll.YData = field(indRec, :);

sound(fieldOnlyNoise(indRec,:), 44100)
sound(field(indRec,:), 44100)


