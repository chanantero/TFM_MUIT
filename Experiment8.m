%% Experiment 8

% Similar to Experiment7.m, but using mainly time processing, even if it's
% slower than the frequency processing.

% The analysis won't be so exhaustive, since the purpose is just proving
% that the time processing works and that the analysis processing would be
% a valid option for analysis because the correspondence between the two
% has been proven.

%% Preamble
pathSetUp;

imagesPath = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\TFM\Img\';

dataPathName = [globalPath, 'Data\'];
ID = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

%% System set up.
% obj = SimulationController;

% Constants
c = 340; % Sound velocity (m/s)
fs = 44100; % Sample frequency (samples/s)
WFSarrayOffset = [0.5, 1.5, 1.5]; % [x, y, z] coordinates. Useful for generating acoustic path IR.

% Noise source coefficient
amplitude = 1;
phase = 0;

% Filter variables for the time WFS filter.
WFSfilterLength = 22050;
% % Creation of frequency filters with different orders.
% magnFiltOrder = 2.^(12);
% hilbertFiltOrder = 2.^(12);
% numFreqFilters = length(magnFiltOrder);
% 
% freqFilters = cell(numFreqFilters, 1);
% freqFiltDelays = zeros(numFreqFilters, 1);
% for k = 1:numFreqFilters
%     [freqFilter, delay] = getFrequencyFilter( magnFiltOrder(k), hilbertFiltOrder(k), fs );    
%     freqFilters{k} = freqFilter;
%     freqFiltDelays(k) = delay;
% end

% Microphone positions
% Rectangular grid
marginRatio = 0.3;
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
recPositions = [X(:), Y(:), Z(:)];
% Comment later
recPositions = [centerX, centerY, 0];

% Positions of the noise source
% Quarter of a circle
numPointsPerQuarter = 4;
radius = [5 10];
numCircles = numel(radius);
alpha = linspace(0, pi/2, numPointsPerQuarter)';
xOctagon = obj.WFSposition(:, 1);
yOctagon = obj.WFSposition(:, 2);
centreX = (max(xOctagon) + min(xOctagon))/2;
centreY = (max(yOctagon) + min(yOctagon))/2;
x = centreX + repmat(radius, numPointsPerQuarter, 1).*repmat(cos(alpha), 1, numCircles);
y = centreY + repmat(radius, numPointsPerQuarter, 1).*repmat(sin(alpha), 1, numCircles);
NSpositions = [x(:), y(:), zeros(numel(x), 1)];
% Comment later
NSpositions = [centreX + 5, centerY, 0];

% Frequencies
freqs = [440]; numFreqs = length(freqs);

% Room characteristics and impulse response of chamber
numReverbTime = 2;
beta = linspace(0, 1, numReverbTime); % Average reflection coefficient of the walls of the chamber
WFS_AcPath_previously_calculated = false;
NS_AcPath_previously_calculated = true;
appendFreeSpaceAcPaths = false;
% Comment later
beta = 0;

% WFS options
frequencyCorrection = false;

% Simulation options
timeDomainActive = true;
fakeTimeProcessing = false;
frequencyDomainActive = true;
automaticLengthModification = false;
predefSignals = true;
saveSignals = true;

durSign = 1; % Duration of tone for time processing
t = (0:ceil(durSign*obj.Fs)-1)/obj.Fs;
NSsignal = chirp(t, 20, durSign, 940);

%% Setup first parameters
SetupParametersScript

%% Pre-calculate impulse responses
% Since it would take a lot of time to calculate the IR during the
% simulations, we calculate it previously and save it to a .mat. 

AcousticPathCalculationScript
% WFSposition = obj.WFSposition;
% save([dataPathName, 'AcPathResponsesWFS_', ID, '.mat'], 'r', 'roomDim', 'numSamp', 'WFSposition', 'freqs', 'WFS_IR', 'WFS_FR')
% save([dataPathName, 'AcPathResponsesNS_', ID, '.mat'], 'r', 'roomDim', 'numSamp', 'NSpositions', 'NS_IR', 'freqs', 'NS_FR')

%% Simulation
simulationScript

%% Save info
rec_signals_NFR = rec_signals; % NFR: No Frequency Correction

%% Now with frequency correction
frequencyCorrection = true;
SetupParametersScript
simulationScript

%% Save info
rec_signals_FR = rec_signals; % FR: Frequency Correction

%% Analysis of results

% Perform FFT
recNS = fft(recNS_signals'); %/fs;

    % No Frequency Correction
recWFS_signals_NFR = rec_signals_NFR - recNS_signals;
recWFS_NFR = fft(recWFS_signals_NFR'); %/fs;
rec_NFR = fft(rec_signals_NFR'); %/fs;
corrFact_NFR = -recNS./recWFS_NFR;

    % Frequency Correction
recWFS_signals_FR = rec_signals_FR - recNS_signals;
recWFS_FR = fft(recWFS_signals_FR'); %/fs;
rec_FR = fft(rec_signals_FR'); %/fs;
corrFact_FR = -recNS./recWFS_FR;

% Visualize
numSamp = size(recNS_signals, 2);
t = (0:numSamp-1)/fs;
f = (0:numSamp-1)/(numSamp/fs);

axFR = axes(figure);
plot(axFR, t, recNS_signals(1,:), t, recWFS_signals_FR(1,:), t, rec_signals_FR(1,:))
axFR.XLabel.String = 'Time (s)';
axFR.YLabel.String = 'Signal (arbitrary units)';
axFR.Title.String = 'With Frequency Filter';
legend(axFR, 'NS', 'WFS', 'Total');

axNFR = axes(figure);
plot(axNFR, t, recNS_signals(1,:), t, recWFS_signals_NFR(1,:), t, rec_signals_NFR(1,:))
axNFR.XLabel.String = 'Time (s)';
axNFR.YLabel.String = 'Signal (arbitrary units)';
axNFR.Title.String = 'Without Frequency Filter';
legend(axNFR, 'NS', 'WFS', 'Total');

axFR_FFT = axes(figure);
plot(axFR_FFT, f, abs(corrFact_FR))
axFR_FFT.XLim = [0, 1000];
axFR_FFT.YLim = [0 4];
axFR_FFT.YLabel.String = '|\Psi|';
yyaxis(axFR_FFT, 'right') % Select with ax.YAxisLocation
plot(axFR_FFT, f, rad2deg(angle(corrFact_FR)));
axFR_FFT.XLabel.String = 'Frequency (Hz)';
axFR_FFT.YLabel.String = 'angle(\Psi) (º)';
axFR_FFT.Title.String = 'Correction Factor \Psi with Frequency Filter';

axNFR_FFT = axes(figure);
plot(axNFR_FFT, f, abs(corrFact_NFR))
axNFR_FFT.XLim = [0, 1000];
axNFR_FFT.YLim = [0 4];
axNFR_FFT.YLabel.String = '|\Psi|';
yyaxis(axNFR_FFT, 'right') % Select with ax.YAxisLocation
plot(axNFR_FFT, f, rad2deg(angle(corrFact_NFR)));
axNFR_FFT.XLabel.String = 'Frequency (Hz)';
axNFR_FFT.YLabel.String = 'angle(\Psi) (º)';
axNFR_FFT.Title.String = 'Correction Factor \Psi without Frequency Filter';

% Visualize frequency filter


%% SVG scenario
viewBox = [-WFSarrayOffset(1) -WFSarrayOffset(2) roomDim(1) roomDim(2)];
NSangles = atan2d(centreY - NSpositions(:,2), centreX - NSpositions(:,1));

objSVG = SVGdrawer('viewBox', viewBox, 'NSpositions', NSpositions,...
    'NSangles', NSangles, 'microSymbol', 'dot', 'microSize', 0.05,...
    'microPositions', recPositions);

name = 'Experiment8_scheme';
objSVG.drawSVG([imagesPath, name, '.svg']);

currentFolder = pwd;
cd(imagesPath); % Needed for inkscape to link svg files properly
system(['inkscape -z "', imagesPath, name, '.svg" --export-pdf="', imagesPath, name, '.pdf"'])
cd(currentFolder)