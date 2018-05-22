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
NSpositions = [centreX + 5, centreY, 0];

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
rec_signals_NFC = rec_signals; % NFR: No Frequency Correction

%% Now with frequency correction
frequencyCorrection = true;
SetupParametersScript
simulationScript

%% Save info
rec_signals_FC = rec_signals; % FR: Frequency Correction

%% Analysis of results

% Perform FFT
NS = fft(x');
recNS = fft(recNS_signals'); %/fs;
H_NS = recNS./NS;

    % No Frequency Correction
recWFS_signals_NFC = rec_signals_NFC - recNS_signals;
recWFS_NFC = fft(recWFS_signals_NFC'); %/fs;
rec_NFC = fft(rec_signals_NFC'); %/fs;
H_WFS_NFC = recWFS_NFC./NS;
corrFact_NFC = -recNS./recWFS_NFC;

    % Frequency Correction
recWFS_signals_FC = rec_signals_FC - recNS_signals;
recWFS_FC = fft(recWFS_signals_FC'); %/fs;
rec_FC = fft(rec_signals_FC'); %/fs;
H_WFS_FC = recWFS_FC./NS;
corrFact_FC = -recNS./recWFS_FC;

% Visualize
numSamp = size(recNS_signals, 2);
t = (0:numSamp-1)/fs;
f = (0:numSamp-1)/(numSamp/fs);

axNS = axes(figure);
plot(axNS, t, x)
axNS.XLabel.String = 'Time (s)';
axNS.YLabel.String = 'Signal (arbitrary units)';
axNS.Title.String = 'Transmitted signal';

% printfig(axNS.Parent, imagesPath, 'Experiment8_NStime', 'eps');

axNS_FFT = axes(figure);
plot(axNS_FFT, f, abs(NS));
axNS_FFT.YLabel.String = 'Magnitude (arbitrary)';
axNS_FFT.XLim = [0, 1000];
yyaxis(axNS_FFT, 'right')
plot(axNS_FFT, f, rad2deg(unwrap(angle(NS))));
axNS_FFT.YLabel.String = 'Phase (º)';
axNS_FFT.XLabel.String = 'Frequency (Hz)';
axNS_FFT.Title.String = 'Spectrum of Transmitted Signal';

% printfig(axNS_FFT.Parent, imagesPath, 'Experiment8_NSfreq', 'eps');

axNFC = axes(figure);
plot(axNFC, t, recNS_signals(1,:), t, recWFS_signals_NFC(1,:), t, rec_signals_NFC(1,:))
axNFC.XLabel.String = 'Time (s)';
axNFC.YLabel.String = 'Signal (arbitrary units)';
axNFC.Title.String = 'Without Frequency Filter';
legend(axNFC, 'NS', 'WFS', 'Total');

% printfig(axNFC.Parent, imagesPath, 'Experiment8_signalTime_NFC', 'eps');

axSignNFC_FFT = axes(figure);
aux = [recNS, recWFS_NFC];
plot(axSignNFC_FFT, f, abs(aux))
axSignNFC_FFT.XLim = [0, 1000];
axSignNFC_FFT.YLabel.String = 'Magnitude (arbitrary units)';
yyaxis(axSignNFC_FFT, 'right') % Select with ax.YAxisLocation
plot(axSignNFC_FFT, f, rad2deg(unwrap(angle(aux))));
axSignNFC_FFT.XLabel.String = 'Frequency (Hz)';
axSignNFC_FFT.YLabel.String = 'Phase (º)';
axSignNFC_FFT.Title.String = 'Spectrum of Received Signals without Frequency Filter';
legend(axSignNFC_FFT, 'NS', 'WFS')

% printfig(axSignNFC_FFT.Parent, imagesPath, 'Experiment8_signalFreq_NFC', 'eps');

axH_NFC_FFT = axes(figure);
aux = [H_NS, H_WFS_NFC];
plot(axH_NFC_FFT, f, abs(aux))
axH_NFC_FFT.XLim = [0, 1000];
axH_NFC_FFT.YLabel.String = 'Magnitude (arbitrary units)';
yyaxis(axH_NFC_FFT, 'right') % Select with ax.YAxisLocation
plot(axH_NFC_FFT, f, rad2deg(unwrap(angle(aux))));
axH_NFC_FFT.XLabel.String = 'Frequency (Hz)';
axH_NFC_FFT.YLabel.String = 'Phase (º)';
axH_NFC_FFT.Title.String = 'Transference Functions without Frequency Filter';
legend(axH_NFC_FFT, 'NS', 'WFS')

% printfig(axH_NFC_FFT.Parent, imagesPath, 'Experiment8_transFunc_NFC', 'eps');

axCorr_NFC_FFT = axes(figure);
plot(axCorr_NFC_FFT, f, abs(corrFact_NFC))
axCorr_NFC_FFT.XLim = [0, 1000];
axCorr_NFC_FFT.YLim = [0 4];
axCorr_NFC_FFT.YLabel.String = '|\Psi|';
yyaxis(axCorr_NFC_FFT, 'right') % Select with ax.YAxisLocation
plot(axCorr_NFC_FFT, f, rad2deg(angle(corrFact_NFC)));
axCorr_NFC_FFT.XLabel.String = 'Frequency (Hz)';
axCorr_NFC_FFT.YLabel.String = 'angle(\Psi) (º)';
axCorr_NFC_FFT.Title.String = 'Correction Factor \Psi without Frequency Filter';

% printfig(axCorr_NFC_FFT.Parent, imagesPath, 'Experiment8_CorrFilt_NFC', 'eps');

axFC = axes(figure);
plot(axFC, t, recNS_signals(1,:), t, recWFS_signals_FC(1,:), t, rec_signals_FC(1,:))
axFC.XLabel.String = 'Time (s)';
axFC.YLabel.String = 'Signal (arbitrary units)';
axFC.Title.String = 'With Frequency Filter';
legend(axFC, 'NS', 'WFS', 'Total');

% printfig(axFC.Parent, imagesPath, 'Experiment8_signalTime_FC', 'eps');

axSignFC_FFT = axes(figure);
aux = [recNS, recWFS_FC];
plot(axSignFC_FFT, f, abs(aux))
axSignFC_FFT.XLim = [0, 1000];
axSignFC_FFT.YLabel.String = 'Magnitude (arbitrary units)';
yyaxis(axSignFC_FFT, 'right') % Select with ax.YAxisLocation
plot(axSignFC_FFT, f, rad2deg(unwrap(angle(aux))));
axSignFC_FFT.XLabel.String = 'Frequency (Hz)';
axSignFC_FFT.YLabel.String = 'Phase (º)';
axSignFC_FFT.Title.String = 'Spectrum of Received Signals with Frequency Filter';
legend(axSignFC_FFT, 'NS', 'WFS')

% printfig(axSignFC_FFT.Parent, imagesPath, 'Experiment8_signalFreq_FC', 'eps');

axH_FC_FFT = axes(figure);
aux = [H_NS, H_WFS_FC];
plot(axH_FC_FFT, f, abs(aux))
axH_FC_FFT.XLim = [0, 1000];
axH_FC_FFT.YLabel.String = 'Magnitude (arbitrary units)';
yyaxis(axH_FC_FFT, 'right') % Select with ax.YAxisLocation
plot(axH_FC_FFT, f, rad2deg(unwrap(angle(aux))));
axH_FC_FFT.XLabel.String = 'Frequency (Hz)';
axH_FC_FFT.YLabel.String = 'Phase (º)';
axH_FC_FFT.Title.String = 'Transference Functions with Frequency Filter';

% printfig(axH_FC_FFT.Parent, imagesPath, 'Experiment8_transFunc_FC', 'eps');

axCorrFC_FFT = axes(figure);
plot(axCorrFC_FFT, f, abs(corrFact_FC))
axCorrFC_FFT.XLim = [0, 1000];
axCorrFC_FFT.YLim = [0 4];
axCorrFC_FFT.YLabel.String = '|\Psi|';
yyaxis(axCorrFC_FFT, 'right') % Select with ax.YAxisLocation
plot(axCorrFC_FFT, f, rad2deg(angle(corrFact_FC)));
axCorrFC_FFT.XLabel.String = 'Frequency (Hz)';
axCorrFC_FFT.YLabel.String = 'angle(\Psi) (º)';
axCorrFC_FFT.Title.String = 'Correction Factor \Psi with Frequency Filter';

% printfig(axCorrFC_FFT.Parent, imagesPath, 'Experiment8_CorrFilt_FC', 'eps');

% Visualize frequency filter

%% Diferentes coeficientes de reflexión de la sala
% Room characteristics and impulse response of chamber
numReverbTime = 5;
beta = linspace(0, 0.4, numReverbTime); % Average reflection coefficient of the walls of the chamber
WFS_AcPath_previously_calculated = false;
NS_AcPath_previously_calculated = true;
appendFreeSpaceAcPaths = false;

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

%% Analysis of results

% Visualize
numSamp = size(recNS_signals, 2);
t = (0:numSamp-1)/fs;

axs = gobjects(numReverbTime, 1);
for rt = 1:numReverbTime
    ax = axes(figure);
    axs(rt) = ax;
    plot(ax, t, recNS_signals(1, :, 1, 1, rt), t, rec_signals(1, :, 1, 1, rt));
    ax.XLabel.String = 'Time (s)';
    ax.YLabel.String = 'Signal (arbitrary units)';
%     ax.YLim = [-1, 1];

%     printfig(ax.Parent, imagesPath, ['Experiment8_sign_rev', num2str(rt)], 'eps');
end

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