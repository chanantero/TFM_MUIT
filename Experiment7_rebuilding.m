%% Experiment 7.

% Find the correction factor needed for different frequencies, noise source
% positions, reverberation times, and receiving points.

% The code is very similar to Experiment1.m, but it will be updated with
% the demands of my tutors made on 11th May of 2018. For doing that, I will
% also use code from Experiment6.m.


%% Preamble
pathSetUp;

imagesPath = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\TFM\Img\';

dataPathName = [globalPath, 'Data\'];
ID = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

%% System set up.
obj = SimulationController;

% Constants
c = 340; % Sound velocity (m/s)
fs = 44100; % Sample frequency (samples/s)
WFSfilterLength = 22050;
WFSarrayOffset = [0.5, 1.5, 1.5]; % [x, y, z] coordinates. Useful for generating acoustic path IR.

% Noise source coefficient
amplitude = 1;
phase = 0;

% Filter variables for the time WFS filter. Creation of frequency filters
% with different orders.
magnFiltOrder = 2.^(12);
hilbertFiltOrder = 2.^(12);
numFreqFilters = length(magnFiltOrder);

freqFilters = cell(numFreqFilters, 1);
freqFiltDelays = zeros(numFreqFilters, 1);
for k = 1:numFreqFilters
    [freqFilter, delay] = getFrequencyFilter( magnFiltOrder(k), hilbertFiltOrder(k), fs );    
    freqFilters{k} = freqFilter;
    freqFiltDelays(k) = delay;
end

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
numPointsPerQuarter = 5;
radius = [5 10 15];
numCircles = numel(radius);
alpha = linspace(0, pi/2, numPointsPerQuarter)';
xOctagon = obj.WFSposition(:, 1);
yOctagon = obj.WFSposition(:, 2);
centreX = (max(xOctagon) + min(xOctagon))/2;
centreY = (max(yOctagon) + min(yOctagon))/2;
x = centreX + repmat(radius, numPointsPerQuarter, 1).*repmat(cos(alpha), 1, numCircles);
y = centreY + repmat(radius, numPointsPerQuarter, 1).*repmat(sin(alpha), 1, numCircles);
NSpositions = [x(:), y(:), zeros(numel(x), 1)];

% Frequencies
freqs = 20:40:1000; numFreqs = length(freqs);

% Room characteristics and impulse response of chamber
numReverbTime = 0;
beta = linspace(0, 1, numReverbTime); % Average reflection coefficient of the walls of the chamber
WFS_AcPath_previously_calculated = true;
NS_AcPath_previously_calculated = true;
appendFreeSpaceAcPaths = true;

% WFS options
frequencyCorrection = false;

% Simulation options
timeDomainActive = false;
fakeTimeProcessing = true;
frequencyDomainActive = true;

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

%% Generate information to save
% dimensionOrder = {'domain', 'frequency', 'NSposition', 'ReverbTime'};
% 
% dataStruct = struct();
% dataStruct.domain = {'frequency', 'time'};
% dataStruct.frequency = freqs;
% dataStruct.NSposition = NSpositions;
% dataStruct.ReverbTime = beta;
% dataStruct.dimensionOrder = dimensionOrder;
% dataStruct.data = s;
% dataStruct.recPositions = recPositions;

% save([dataPathName, 'Experiment7data_', ID, '.mat'], 'dataStruct')

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

%% Compare global cancellation for different correction factors
% For each NSposition, you get a global cancellation with the optimum correction 
% factor calculated, and another one with the theoretical one. It is
% obvious that the cancellation with the theoretical correction factor will
% be worse. But, how worse?
freqsMat = pointWiseExtend(freqs, s);
corrFactTheo = sqrt(1i * freqsMat/c);

Cg_theo = zeros(size(s));
Cg_opt = zeros(size(s));
for k = 1:numel(s)
    recCoef_theo = s(k).recNScoef + s(k).recWFScoef * corrFactTheo(k);
    Cg_theo(k) = sum(abs(recCoef_theo).^2)/sum(abs(s(k).recNScoef).^2);
    
    recCoef_opt = s(k).recNScoef + s(k).recWFScoef * corrFactGlobal(k);
    Cg_opt(k) = sum(abs(recCoef_opt).^2)/sum(abs(s(k).recNScoef).^2);
end
Cg_theo_dB = 10*log10(Cg_theo);
Cg_opt_dB = 10*log10(Cg_opt);

% Visualize
freqMat_FS = pointWiseExtend(freqs, s(1, :, :, 1));

Cg_theo_dB_FS = Cg_theo_dB(2, :, :, 1); % FS: Free Space
Cg_opt_dB_FS = Cg_opt_dB(2, :, :, 1); % FS: Free Space

freqStep = freqs(2) - freqs(1);
freqEdges = [(freqs(1:end-1) + freqs(2:end))/2, freqs(end) + freqStep/2];
CancEdges = -20:0;

N_theo = histcounts2(freqMat_FS, Cg_theo_dB_FS, freqEdges, CancEdges);
Nnorm_theo = N_theo./repmat(sum(N_theo, 2), [1, size(N_theo, 2)]);

[N_opt] = histcounts2(freqMat_FS, Cg_opt_dB_FS, freqEdges, CancEdges);
Nnorm_opt = N_opt./repmat(sum(N_opt, 2), [1, size(N_opt, 2)]);

penaltyEdges = 0:15;
[N_diff, freqAux, DifCancEdges] = histcounts2(freqMat_FS, Cg_theo_dB_FS - Cg_opt_dB_FS, freqEdges, penaltyEdges);
Nnorm_diff = N_diff./repmat(sum(N_diff, 2), [1, size(N_diff, 2)]);

ax1 = axes(figure);
C = zeros(length(freqEdges), length(CancEdges));
C(1:end-1, 1:end-1) = Nnorm_theo;
pcolor(ax1, freqEdges, CancEdges, C');
ax1.Title.String = 'C_{global} for theoretical correction factor';
ax1.XLabel.String = 'Frequency (Hz)';
ax1.YLabel.String = 'C_{global}';
colorbar(ax1)
ax1.CLim = [0, 0.5];
% colormap(ax1, 'gray')
% printfig(ax1.Parent, imagesPath, 'Experiment7_CancGlobCorrFactTheo', 'eps');

ax2 = axes(figure);
C = zeros(length(freqEdges), length(CancEdges));
C(1:end-1, 1:end-1) = Nnorm_opt;
p = pcolor(ax2, freqEdges, CancEdges, C');
ax2.Title.String = 'C_{global} for optimum correction factor';
ax2.XLabel.String = 'Frequency (Hz)';
ax2.YLabel.String = 'C_{global}';
colorbar(ax2)
ax2.CLim = [0, 0.5];
% colormap(ax2, 'gray')
% printfig(ax2.Parent, imagesPath, 'Experiment7_CancGlobCorrFactOpt', 'eps');

ax3 = axes(figure);
C = zeros(length(freqEdges), length(DifCancEdges));
C(1:end-1, 1:end-1) = Nnorm_diff;
p = pcolor(ax3, freqEdges, DifCancEdges, C');
ax3.Title.String = 'Cancellation penalty';
ax3.XLabel.String = 'Frequency (Hz)';
ax3.YLabel.String = '$10\log_{10}(C_{global,theo}) - 10\log_{10}(C_{global,opt})$';
ax3.YLabel.Interpreter = 'latex';
colorbar(ax3)
% colormap(ax3, 'gray')
% printfig(ax3.Parent, imagesPath, 'Experiment7_CancGlobPenalty', 'eps');

%% Increase complexity of correction factor step by step

% A) Scalation by a unique real coefficient.
% B) Scalation by a unique complex coefficient
refFreq = 440;
corrFact_A = sqrt(refFreq/c);
corrFact_B = sqrt(1i*refFreq/c);

freqsMat = pointWiseExtend(freqs, s);

Cg_A = zeros(size(s));
Cg_B = zeros(size(s));

for k = 1:numel(s)
    recCoef_A = s(k).recWFScoef * corrFact_A + s(k).recNScoef;
    recCoef_B = s(k).recWFScoef * corrFact_B + s(k).recNScoef;

    Cg_A(k) = sum(abs(recCoef_A).^2)/sum(abs(s(k).recNScoef).^2);
    Cg_B(k) = sum(abs(recCoef_B).^2)/sum(abs(s(k).recNScoef).^2);
end
Cg_A_dB = 10*log10(Cg_A);
Cg_B_dB = 10*log10(Cg_B);

% Visualize
freqMat_FS = pointWiseExtend(freqs, s(1, :, :, 1));

Cg_A_dB_FS = Cg_A_dB(2, :, :, 1); % FS: Free Space
Cg_B_dB_FS = Cg_B_dB(2, :, :, 1); % FS: Free Space

freqStep = freqs(2) - freqs(1);
freqEdges = [0, (freqs(1:end-1) + freqs(2:end))/2, freqs(end) + freqStep/2];
CancEdges = -20:0;

N_A = histcounts2(freqMat_FS, Cg_A_dB_FS, freqEdges, CancEdges);
Nnorm_A = N_A./repmat(sum(N_A, 2), [1, size(N_A, 2)]); Nnorm_A(isnan(Nnorm_A)) = 0;

N_B = histcounts2(freqMat_FS, Cg_B_dB_FS, freqEdges, CancEdges);
Nnorm_B = N_B./repmat(sum(N_B, 2), [1, size(N_B, 2)]); Nnorm_B(isnan(Nnorm_B)) = 0;

ax1 = axes(figure);
C = zeros(length(freqEdges), length(CancEdges));
C(1:end-1, 1:end-1) = Nnorm_A;
pcolor(ax1, freqEdges, CancEdges, C');
ax1.Title.String = ['$C_{global}$ for $\Psi=\sqrt{\frac{', num2str(refFreq), '}{c}}$'];
ax1.Title.Interpreter = 'latex';
ax1.XLabel.String = 'Frequency (Hz)';
ax1.YLabel.String = '$C_{global}$';
ax1.YLabel.Interpreter = 'latex';
colorbar(ax1)
% printfig(ax1.Parent, imagesPath, 'Experiment7_CancGlobCorrFactA', 'eps');

ax2 = axes(figure);
C = zeros(length(freqEdges), length(CancEdges));
C(1:end-1, 1:end-1) = Nnorm_B;
pcolor(ax2, freqEdges, CancEdges, C');
ax2.Title.String = ['$C_{global}$ for $\Psi=\sqrt{j\frac{', num2str(refFreq), '}{c}}$'];
ax2.Title.Interpreter = 'latex';
ax2.XLabel.String = 'Frequency (Hz)';
ax2.YLabel.String = '$C_{global}$';
ax2.YLabel.Interpreter = 'latex';
colorbar(ax2)
% printfig(ax2.Parent, imagesPath, 'Experiment7_CancGlobCorrFactB', 'eps');

%% Continuation of previous section: simulation of full correction factor in time domain

% Creation of frequency filters with different orders
magnFiltOrder = 2.^(9:12);
hilbertFiltOrder = 2.^(9:12);
numFreqFilters = length(magnFiltOrder);

freqFilters = cell(numFreqFilters, 1);
freqFiltDelays = zeros(numFreqFilters, 1);
for k = 1:numFreqFilters
    [freqFilter, delay] = getFrequencyFilter( magnFiltOrder(k), hilbertFiltOrder(k), fs );    
    freqFilters{k} = freqFilter;
    freqFiltDelays(k) = delay;
end

% WFS options
frequencyCorrection = true;

% Simulation options
timeDomainActive = true;
fakeTimeProcessing = true;
frequencyDomainActive = true;

%% Setup first parameters
SetupParametersScript

%% Simulation
simulationScript;

%% Global cancellation with the filter response

Cg = zeros(size(s));
for k = 1:numel(s)
    Cg(k) = sum(abs(s(k).recCoef).^2)/sum(abs(s(k).recNScoef).^2);
end
Cg_dB = 10*log10(Cg);

% Visualize
freqMat_FS = pointWiseExtend(freqs, s(1, :, :, 1));

freqStep = freqs(2) - freqs(1);
freqEdges = [(freqs(1:end-1) + freqs(2:end))/2, freqs(end) + freqStep/2];
CancEdges = -20:0;

axs = gobjects(numFreqFilters + 1, 1);
for filt = 1:numFreqFilters + 1
    Cg_dB_FS = Cg_dB(filt, :, :, 1); % FS: Free Space
    N = histcounts2(freqMat_FS, Cg_dB_FS, freqEdges, CancEdges);
    Nnorm = N./repmat(sum(N, 2), [1, size(N, 2)]); Nnorm(isnan(Nnorm)) = 0;
    
    ax = axes(figure);
    axs(filt) = ax;
    C = zeros(length(freqEdges), length(CancEdges));
    C(1:end-1, 1:end-1) = Nnorm;
    pcolor(ax, freqEdges, CancEdges, C');
    ax.Title.String = 'Global Cancellation';
    ax.Title.Interpreter = 'latex';
    ax.XLabel.String = 'Frequency (Hz)';
    ax.YLabel.String = '$C_{global}$';
    ax.YLabel.Interpreter = 'latex';
    colorbar(ax)
    
    name = ['Experiment7_CancGlobFreqFilt_', num2str(filt)];        
%     printfig(ax.Parent, imagesPath, name, 'eps');
end

%% Different reverberation times, perfect frequency filter
numReverbTime = 5;
beta = linspace(0, 1, numReverbTime);
WFS_AcPath_previously_calculated = false;
NS_AcPath_previously_calculated = true;
appendFreeSpaceAcPaths = false;

% WFS options
frequencyCorrection = false;

% Simulation options
timeDomainActive = false;
frequencyDomainActive = true;

freqFilters = [];

SetupParametersScript
AcousticPathCalculationScript % Pre-calculate impulse responses

%% Simulation
simulationScript;

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

% Visualize
freqs_aux = 0:1000;
correctFactTheo = sqrt(1i * freqs_aux/c);

corrFactGlobal_zeroReverb = corrFactGlobal(end, :, :, 1);
freqMat_glob = pointWiseExtend(freqs, corrFactGlobal_zeroReverb);

freqStep = freqs(2) - freqs(1);
freqEdges = [(freqs(1:end-1) + freqs(2:end))/2, freqs(end) + freqStep/2];
magCorrFactEdges = 0:0.1:3;
phaseCorrFactEdges = -180:4:180;

% [N_abs_glob, Xedges_abs_glob, Yedges_abs_glob] = ...
%     histcounts2(freqMat_glob, abs(corrFactGlobal_zeroReverb), freqEdges, magCorrFactEdges);
% [N_phase_glob, Xedges_phase_glob, Yedges_phase_glob] = ...
%     histcounts2(freqMat_glob, rad2deg(angle(corrFactGlobal_zeroReverb)), freqEdges, magCorrFactEdges);
% 
% Nnorm_abs_glob = N_abs_glob./repmat(sum(N_abs_glob, 2), [1, size(N_abs_glob, 2)]);
% Nnorm_phase_glob = N_phase_glob./repmat(sum(N_phase_glob, 2), [1, size(N_phase_glob, 2)]);

% ax1 = axes(figure);
% C = zeros(length(Xedges_abs_glob), length(Yedges_abs_glob));
% C(1:end-1, 1:end-1) = Nnorm_abs_glob;
% pcolor(ax1, Xedges_abs_glob, Yedges_abs_glob, C');
% ax1.Title.String = '|\Psi_{global}| for global cancellation';
% ax1.XLabel.String = 'Frequency (Hz)';
% ax1.YLabel.String = '|\Psi_{global}|';
% colorbar
% ax1.NextPlot = 'Add';
% plot(ax1, freqs_aux, abs(correctFactTheo), 'r', 'LineWidth', 4);
% 
% % printfig(ax3.Parent, imagesPath, 'Experiment7_corrFactGlobCancAbs', 'eps');
% 
% ax2 = axes(figure);
% C = zeros(length(Xedges_phase_glob), length(Yedges_phase_glob));
% C(1:end-1, 1:end-1) = Nnorm_phase_glob;
% pcolor(ax2, Xedges_phase_glob, Yedges_phase_glob, C');
% ax2.Title.String = 'phase(\Psi_{global}) for global cancellation';
% ax2.XLabel.String = 'Frequency (Hz)';
% ax2.YLabel.String = 'phase(\Psi_{global})';
% colorbar

% printfig(ax4.Parent, imagesPath, 'Experiment7_corrFactGlobCancPhase', 'eps');

axs = gobjects(numReverbTime, 1);
for rt = 1:numReverbTime
    N = histcounts2(freqMat_FS, abs(corrFactGlobal(end, :, :, rt)), freqEdges, magCorrFactEdges);
    Nnorm = N./repmat(sum(N, 2), [1, size(N, 2)]); Nnorm(isnan(Nnorm)) = 0;
    
    ax = axes(figure);
    axs(rt) = ax;
    C = zeros(length(freqEdges), length(magCorrFactEdges));
    C(1:end-1, 1:end-1) = Nnorm;
    pcolor(ax, freqEdges, magCorrFactEdges, C');
    ax.Title.String = '$|\Psi|$';
    ax.Title.Interpreter = 'latex';
    ax.XLabel.String = 'Frequency (Hz)';
    ax.YLabel.String = '$|\Psi|$';
    ax.YLabel.Interpreter = 'latex';
    colorbar(ax)
    ax.NextPlot = 'Add';
    plot(ax, freqs_aux, abs(correctFactTheo), 'r', 'LineWidth', 4);
    
    name = ['Experiment7_corrFactGlob_rev', num2str(rt)];        
%     printfig(ax.Parent, imagesPath, name, 'eps');
end

%% SVG scenario
viewBox = [-WFSarrayOffset(1) -WFSarrayOffset(2) roomDim(1) roomDim(2)];
NSangles = atan2d(centreY - NSpositions(:,2), centreX - NSpositions(:,1));

objSVG = SVGdrawer('viewBox', viewBox, 'NSpositions', NSpositions,...
    'NSangles', NSangles, 'microSymbol', 'dot', 'microSize', 0.05,...
    'microPositions', recPositions);

name = 'Experiment7_scheme';
objSVG.drawSVG([imagesPath, name, '.svg']);

currentFolder = pwd;
cd(imagesPath); % Needed for inkscape to link svg files properly
system(['inkscape -z "', imagesPath, name, '.svg" --export-pdf="', imagesPath, name, '.pdf"'])
cd(currentFolder)


