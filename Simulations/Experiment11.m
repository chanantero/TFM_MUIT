%% Experiment 11.

% Based on Experiment7_rebuilding.m.
% Created on 27th June 2018, it aims at generating the graphs for the
% introduction to simulations chapter in the thesis report.

%% Preamble
pathSetUp;

imagesPath = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\TFM\Img\';

dataPathName = [globalPath, 'Data\'];
ID = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

%% System set up.
if ~exist('obj', 'var') || ~isvalid(obj)
    obj = SimulationController;
end

% Constants
c = 340; % Sound velocity (m/s)
fs = 44100; % Sample frequency (samples/s)
WFSfilterLength = 22050;
WFSarrayOffset = [0.5, 1.5, 1.65]; % [x, y, z] coordinates. Useful for generating acoustic path IR.

% Noise source coefficient
amplitude = 1;
phase = 0;

% Filter variables for the time WFS filter. Creation of frequency filters
% with different orders.
magnFiltOrder = 2.^(12);
hilbertFiltOrder = 2.^(12);
[freqFilter, freqFiltDelay] = getFrequencyFilter( magnFiltOrder, hilbertFiltOrder, fs );    
freqFilters = {freqFilter};
freqFiltDelays = freqFiltDelay;

% Microphone positions
% Rectangular grid
marginRatio = 0.6;
numPointsX = 5;
numPoinstY = 5;
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

% Positions of the noise source
% Quarter of a circle
NSpositions = [centerX + 5, centerY, 0];

% Frequencies
freqs = 440;

% Room characteristics and impulse response of chamber
beta = 0; % Average reflection coefficient of the walls of the chamber
WFS_AcPath_previously_calculated = true;
NS_AcPath_previously_calculated = true;
appendFreeSpaceAcPaths = true;

% WFS options
frequencyCorrection = false;
attenuationType = 'Ruben';

% Simulation options
timeDomainActive = true;
fakeTimeProcessing = false;
frequencyDomainActive = true;
automaticLengthModification = false;
saveSignals = true;

% Duration of tone for time processing
durSign = 1;

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

%% Graphs for the report

% Time domain
numSamples = size(rec_signals, 2);
t = (0:numSamples - 1)/fs;
ax = axes(figure);
plot(ax, t, recNS_signals(13, :), t, rec_signals(13, :))
ax.XLim = [0.2, 0.21];
ax.XLabel.String = 'Time (s)';
printfig(ax.Parent, imagesPath, 'Experiment11_exampleNoFreqCorr', 'eps');



%% Load saved information
% s = dataStruct.data;
% beta = dataStruct.ReverbTime;
% 
% numMicro = size(dataStruct.recPositions, 1);


%% Visualization: 2D map, case by case

% Create simulationViewer object
objVis = simulationViewer(obj.ax, s);

%% Visualization: correction factor
% Calculate individual and global correction factors
corrFactInd = zeros([size(s), numMicro]); % Individual correction factor
numDims = ndims(s);
subinds_s = cell(numDims, 1);
subinds_ext = cell(numDims + 1, 1);
for k = 1:numel(s)
    s(k).corrFactIndividual = -s(k).recNScoef./s(k).recWFScoef;
    s(k).corrFactGlobal = -s(k).recWFScoef\s(k).recNScoef;
    
    [subinds_s{:}] = ind2sub(size(s), k);
    for d = 1:numDims
        subinds_ext{d} = subinds_s{d}*ones(numMicro, 1);
    end
    subinds_ext{end} = (1:numMicro)';
    inds = sub2ind(size(corrFactInd), subinds_ext{:})';
    
    corrFactInd(inds) = s(k).corrFactIndividual;
end
corrFactGlobal = zeros(size(s)); % Individual correction factor
corrFactGlobal(:) = [s.corrFactGlobal];

% How disperse are the individual correction factors for different NS
% positions and microphone positions?
% Note: you can use histogram2 instead of pcolor for visualizing, but the
% flexibility is smaller
freqs_aux = 0:1000;
correctFactTheo = sqrt(1i * freqs_aux/c);

corrFactInd_zeroReverb = corrFactInd(2, :, :, 1, :);
corrFactGlobal_zeroReverb = corrFactGlobal(2, :, :, 1);
freqMat_ind = pointWiseExtend(freqs, corrFactInd_zeroReverb);
freqMat_glob = pointWiseExtend(freqs, corrFactGlobal_zeroReverb);
[N_abs, Xedges_abs, Yedges_abs] = histcounts2(freqMat_ind, abs(corrFactInd_zeroReverb), 'BinWidth', [freqs(2) - freqs(1), 0.1]);
[N_phase, Xedges_phase, Yedges_phase] = histcounts2(freqMat_ind, rad2deg(angle(corrFactInd_zeroReverb)), 'BinWidth', [freqs(2) - freqs(1), 1]);
[N_abs_glob, Xedges_abs_glob, Yedges_abs_glob] = ...
    histcounts2(freqMat_glob, abs(corrFactGlobal_zeroReverb), 'BinWidth', [freqs(2) - freqs(1), 0.1]);
[N_phase_glob, Xedges_phase_glob, Yedges_phase_glob] = ...
    histcounts2(freqMat_glob, rad2deg(angle(corrFactGlobal_zeroReverb)), 'BinWidth', [freqs(2) - freqs(1), 1]);

Nnorm_abs = N_abs./repmat(sum(N_abs, 2), [1, size(N_abs, 2)]);
Nnorm_phase = N_phase./repmat(sum(N_phase, 2), [1, size(N_phase, 2)]);
Nnorm_abs_glob = N_abs_glob./repmat(sum(N_abs_glob, 2), [1, size(N_abs_glob, 2)]);
Nnorm_phase_glob = N_phase_glob./repmat(sum(N_phase_glob, 2), [1, size(N_phase_glob, 2)]);

ax1 = axes(figure);
C = zeros(length(Xedges_abs), length(Yedges_abs));
C(1:end-1, 1:end-1) = Nnorm_abs;
pcolor(ax1, Xedges_abs, Yedges_abs, C');
ax1.Title.String = '|\Psi_{ind}| for individual cancellation';
ax1.XLabel.String = 'Frequency (Hz)';
ax1.YLabel.String = '|\Psi_{ind}|';
colorbar
ax1.NextPlot = 'Add';
plot(ax1, freqs_aux, abs(correctFactTheo), 'r', 'LineWidth', 4);

% printfig(ax1.Parent, imagesPath, 'Experiment7_corrFactIndCancAbs', 'eps');

ax2 = axes(figure);
C = zeros(length(Xedges_phase), length(Yedges_phase));
C(1:end-1, 1:end-1) = Nnorm_phase;
pcolor(ax2, Xedges_phase, Yedges_phase, C');
ax2.Title.String = 'phase(\Psi_{ind}) for individual cancellation';
ax2.XLabel.String = 'Frequency (Hz)';
ax2.YLabel.String = 'phase(\Psi_{ind})';
colorbar

% printfig(ax2.Parent, imagesPath, 'Experiment7_corrFactIndCancPhase', 'eps');

ax3 = axes(figure);
C = zeros(length(Xedges_abs_glob), length(Yedges_abs_glob));
C(1:end-1, 1:end-1) = Nnorm_abs_glob;
pcolor(ax3, Xedges_abs_glob, Yedges_abs_glob, C');
ax3.Title.String = '|\Psi_{global}| for global cancellation';
ax3.XLabel.String = 'Frequency (Hz)';
ax3.YLabel.String = '|\Psi_{global}|';
colorbar
ax3.NextPlot = 'Add';
plot(ax3, freqs_aux, abs(correctFactTheo), 'r', 'LineWidth', 4);

% printfig(ax3.Parent, imagesPath, 'Experiment7_corrFactGlobCancAbs', 'eps');

ax4 = axes(figure);
C = zeros(length(Xedges_phase_glob), length(Yedges_phase_glob));
C(1:end-1, 1:end-1) = Nnorm_phase_glob;
pcolor(ax4, Xedges_phase_glob, Yedges_phase_glob, C');
ax4.Title.String = 'phase(\Psi_{global}) for global cancellation';
ax4.XLabel.String = 'Frequency (Hz)';
ax4.YLabel.String = 'phase(\Psi_{global})';
colorbar

% printfig(ax4.Parent, imagesPath, 'Experiment7_corrFactGlobCancPhase', 'eps');
