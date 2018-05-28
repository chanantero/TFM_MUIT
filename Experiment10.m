%% Experiment 10
% Comparison of the WFS attenuations used by Miguel and by me.

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
d = 0.18; % Separation between WFS array loudspeakers

% Noise source coefficient
amplitude = 1;
phase = 0;

% Filter variables for the time WFS filter.
freqFilters = {1}; % Required for technical reasons
freqFiltDelays = 0; % Required for technical reasons

% Frequencies
freqs = 0; % Random frequency. Required for technical purposes
numFreqs = 1; % Required for technical reasons

% Room characteristics and impulse response of chamber
WFS_AcPath_previously_calculated = false;
NS_AcPath_previously_calculated = true;
appendFreeSpaceAcPaths = false;
numSampIR = 1024;
predefNumSampIR = true;
numReverbTime = 1;
beta = 0;

% WFS options
frequencyCorrection = false;
attenuationType = 'Miguel';
roomDim = [9.1300    4.4800    2.6400];

% Simulation options
timeDomainActive = true;
fakeTimeProcessing = false;
frequencyDomainActive = false;
automaticLengthModification = false;
predefSignals = true;
saveSignals = true;

zPos = 1.65;
NSpositions = [0.1 2.3 zPos];

predefRoomDim = true;
[X, Y] = ndgrid(mallado_x, mallado_y);
recPositions = [X(:), Y(:), zPos*ones(numel(X), 1)];

WFSposition = [WFSpos', zPos*ones(96, 1)];
obj.WFSposition = WFSposition;
broadsideDir = [cosd(tecta' - 90), sind(tecta' - 90), zeros(96, 1)];
obj.WFSToolObj.WFSarrayOrientation = simulator.vec2rotVec(broadsideDir);
obj.WFSToolObj.scenarioObj.ax.XLim = [0, roomDim(1)];
obj.WFSToolObj.scenarioObj.ax.YLim = [0, roomDim(2)];

WFSfilterLength = size(filtros_array, 1);
predefWFSfilterLength = true;

WFSarrayOffset = [0, 0, 0]; % [x, y, z] coordinates. Useful for generating acoustic path IR.

SetupParametersScript;
AcousticPathCalculationScript;

simulationScript;

%% Analysis of results

% Visualize
    % My code
recWFS_signals_RubScript = rec_signals - recNS_signals;

    % Miguel's code
recWFS_signals_Miguel = [POT_ar', zeros(numMicro, length(x) - size(POT_ar, 1))]; % Mucho cuidado con esto!!! Al principio has interpretado que POT_ar era el campo recibido, pero es el campo proveniente del WFS!! Esto es importante. Te ha hecho perder varias horas
recNS_signals_Miguel = [POT_ad', zeros(numMicro, length(x) - size(POT_ad, 1))];

% Visualize
numSamp = size(recNS_signals, 2);
t = (0:numSamp-1)/fs;

axNS = axes(figure);
plot(axNS, t, recNS_signals(1,:), t, recNS_signals_Miguel(1,:))
axNS.XLabel.String = 'Time (s)';
axNS.YLabel.String = 'Signal (arbitrary units)';
axNS.Title.String = 'Received signal from the noise source';
plot(axNS, t, recNS_signals')

axWFS = axes(figure);
plot(axWFS, t, recWFS_signals(1,:), t, d*recWFS_signals_Miguel(1,:), t, recNS_signals_Miguel(1,:))
axWFS.XLabel.String = 'Time (s)';
axWFS.YLabel.String = 'Signal (arbitrary units)';
axWFS.Title.String = 'Received signal from the WFS array';

recNScoef = signal2pulseCoefficientMatrix([0 durSign], freqCos, 1, (recNS_signals)', fs);
recWFScoef = signal2pulseCoefficientMatrix([0 durSign], freqCos, 1, (recWFS_signals)', fs);
corrFact = -recNScoef./recWFScoef;
abs(corrFact)
rad2deg(angle(corrFact))

% Why is the received signal from WFS slightly different?
% Answer: the WFS filters are different. Some of your deltas are delayed
% one sample with respect the deltas of Miguel's code. This happens because I use the
% function floor and he uses round to calculate the sample where de delta
% will be located

% Compare WFS filters. You can see differences of one sample in some deltas
tWFSfilt = (0:WFSfilterLength - 1)/fs;
ax = axes(figure);
wfsInd = 2;
plot(ax, tWFSfilt, filtros_array(:, wfsInd), tWFSfilt, permute(obj.WFSToolObj.filtersWFS_IR(wfsInd, 1, :), [3 1 2]))

% Let's test if that's the only remaining difference
obj.WFSToolObj.filtersWFS_IR = repmat(permute(filtros_array, [2 3 1]), [1 2 1]);
isequal(repmat(permute(filtros_array, [2 3 1]), [1 2 1]), obj.WFSToolObj.filtersWFS_IR)
obj.WFSToolObj.WFScalculation();
obj.WFSToolObj.simulate();
rec_signal2 = obj.WFSToolObj.simulField;
recWFS_signals2 = rec_signal2 - recNS_signals;

axWFS = axes(figure);
plot(axWFS, t, recWFS_signals2(1,:), t, recWFS_signals_Miguel(1,:))
axWFS.XLabel.String = 'Time (s)';
axWFS.YLabel.String = 'Signal (arbitrary units)';
axWFS.Title.String = 'Received signal from the WFS array';

% Conclussion: my way of doing it and Miguel's way are the same.
% Hence, it is confirmed that the only difference between my code and his
% code is that he uses another r0 and doesn't multiply by the loudspeaker
% spearation.

%% Analysis of the two attenuation calculations: mine and Miguel's
% One position (the same as previously)
% Using just my scripts, because it has been proved the correspondence
% between Miguel's script and mine
% Same sinusoidal signal at one frequency

attenuationType = 'Ruben';
SetupParametersScript;
simulationScript;
rec_signals_Rub = rec_signals;
recNS_signals_Rub = recNS_signals;

attenuationType = 'Miguel';
SetupParametersScript;
simulationScript;
rec_signals_Mig = rec_signals;
recNS_signals_Mig = recNS_signals;

assert(isequal(recNS_signals_Mig, recNS_signals_Rub), 'Received field from noise source should be the same!')
recNS_signals = recNS_signals_Rub;

% Visualize
    % My code
recWFS_signals_Rub = rec_signals_Rub - recNS_signals;
recWFS_signals_Mig = rec_signals_Mig - recNS_signals;
isequal(recWFS_signals_Mig, recWFS_signals_RubScript)

% Visualize
numSamp = size(recNS_signals, 2);
t = (0:numSamp-1)/fs;

axWFS = axes(figure);
plot(axWFS, t, recWFS_signals_Rub(1,:), t, recWFS_signals_Mig(1,:), t, recNS_signals(1,:))
axWFS.XLabel.String = 'Time (s)';
axWFS.YLabel.String = 'Signal (arbitrary units)';
axWFS.Title.String = 'Received signal from the WFS array';

recNScoef = signal2pulseCoefficientMatrix([0 durSign], freqCos, 1, (recNS_signals)', fs);
recWFScoef_Rub = signal2pulseCoefficientMatrix([0 durSign], freqCos, 1, (recWFS_signals_Rub)', fs);
recWFScoef_Mig = signal2pulseCoefficientMatrix([0 durSign], freqCos, 1, (recWFS_signals_Mig)', fs);
corrFact_Rub = -recNScoef./recWFScoef_Rub;
corrFact_Mig = -recNScoef./recWFScoef_Mig;
abs(corrFact_Rub)
rad2deg(angle(corrFact_Rub))
abs(corrFact_Mig)
rad2deg(angle(corrFact_Mig))

%% Analysis of the two attenuation calculations, over the bandwith
% One position (the same as previously)
% Using just my scripts, because it has been proved the correspondence
% between Miguel's script and mine
% Chirp signal for the noise source
durSign = 1; % Duration of tone for time processing
t = (0:ceil(durSign*obj.Fs)-1)/obj.Fs;
NSsignal = chirp(t, 20, durSign, 940);

WFSfilterLength = 2^12;
numSampIR = 2^12;
WFS_AcPath_previously_calculated = true;
SetupParametersScript;
AcousticPathCalculationScript

attenuationType = 'Ruben';
SetupParametersScript;
simulationScript;
rec_signals_Rub = rec_signals;
recNS_signals_Rub = recNS_signals;

attenuationType = 'Miguel';
SetupParametersScript;
simulationScript;
rec_signals_Mig = rec_signals;
recNS_signals_Mig = recNS_signals;

% Analysis
assert(isequal(recNS_signals_Mig, recNS_signals_Rub), 'Signals from NS should be the same')
recNS_time = recNS_signals_Mig;

% Perform FFT
recNS_freq = fft(permute(recNS_time, [2 1 3 4]));

    % Rubén attenuation
recWFS_time_Rub = rec_signals_Rub - recNS_time;
recWFS_freq_Rub = fft(permute(recWFS_time_Rub, [2 1 3 4])); %/fs;
rec_freq_Rub = fft(permute(rec_signals_Rub, [2 1 3 4])); %/fs;
corrFact_Rub = -recNS_freq./recWFS_freq_Rub;

    % Miguel's attenuation
recWFS_time_Mig = rec_signals_Mig - recNS_time;
recWFS_freq_Mig = fft(permute(recWFS_time_Mig, [2 1 3 4])); %/fs;
rec_freq_Mig = fft(permute(rec_signals_Mig, [2 1 3 4])); %/fs;
corrFact_Mig = -recNS_freq./recWFS_freq_Mig;

numSamp = size(recNS_time, 2);
f = (0:numSamp - 1)*fs/numSamp;
corrFact_Rub_phase = rad2deg(angle(corrFact_Rub));
corrFact_Rub_abs = abs(corrFact_Rub);
ax = axes(figure);
ind = f >=0 & f<1000;
plot(ax, f(ind), corrFact_Rub_phase(ind, 1))

%% Analysis of the two attenuation calculations. Chirp and multiple positions

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
NSpositions = [x(:), y(:), zPos*ones(numel(x), 1)];
NSpositions = [; NSpositions];

durSign = 1; % Duration of tone for time processing
t = (0:ceil(durSign*obj.Fs)-1)/obj.Fs;
NSsignal = chirp(t, 20, durSign, 940);

WFSfilterLength = 2^12;
numSampIR = 2^12;

WFS_AcPath_previously_calculated = true;
attenuationType = 'Miguel';
SetupParametersScript
AcousticPathCalculationScript
simulationScript;
rec_signals_Mig = rec_signals;
recNS_signals_Mig = recNS_signals;

attenuationType = 'Ruben';
SetupParametersScript;
simulationScript;
rec_signals_Rub = rec_signals;
recNS_signals_Rub = recNS_signals;


%% Analysis of results

% Visualize

assert(isequal(recNS_signals_Mig, recNS_signals_Rub), 'Signals from NS should be the same')
recNS_time = recNS_signals_Mig;

% Perform FFT
recNS_freq = fft(permute(recNS_time, [2 1 3 4]));

    % Rubén attenuation
recWFS_time_Rub = rec_signals_Rub - recNS_time;
recWFS_freq_Rub = fft(permute(recWFS_time_Rub, [2 1 3 4])); %/fs;
rec_freq_Rub = fft(permute(rec_signals_Rub, [2 1 3 4])); %/fs;
corrFact_Rub = -recNS_freq./recWFS_freq_Rub;

    % Miguel's attenuation
recWFS_time_Mig = rec_signals_Mig - recNS_time;
recWFS_freq_Mig = fft(permute(recWFS_time_Mig, [2 1 3 4])); %/fs;
rec_freq_Mig = fft(permute(rec_signals_Mig, [2 1 3 4])); %/fs;
corrFact_Mig = -recNS_freq./recWFS_freq_Mig;

% Graphs


% Histograms
numSamp = size(recNS_time, 2);
t = (0:numSamp-1)/fs;
f = (0:numSamp-1)*fs/numSamp;
freqEdges = 10:10:1000;
absCorrFactEdges = 0:0.1:4;
phaseCorrFactEdges = 0:90;

freqMat_ind = pointWiseExtend(f', corrFact_Rub);
N_abs_Rub = histcounts2(freqMat_ind, abs(corrFact_Rub), freqEdges, absCorrFactEdges);
N_phase_Rub = histcounts2(freqMat_ind, rad2deg(angle(corrFact_Rub)), freqEdges, phaseCorrFactEdges);
N_abs_Mig = histcounts2(freqMat_ind, abs(corrFact_Mig), freqEdges, absCorrFactEdges);
N_phase_Mig = histcounts2(freqMat_ind, rad2deg(angle(corrFact_Mig)), freqEdges, phaseCorrFactEdges);

Nnorm_abs_Rub = N_abs_Rub./repmat(sum(N_abs_Rub, 2), [1, size(N_abs_Rub, 2)]);
Nnorm_phase_Rub = N_phase_Rub./repmat(sum(N_phase_Rub, 2), [1, size(N_phase_Rub, 2)]);
Nnorm_abs_Mig = N_abs_Mig./repmat(sum(N_abs_Mig, 2), [1, size(N_abs_Mig, 2)]);
Nnorm_phase_Mig = N_phase_Mig./repmat(sum(N_phase_Mig, 2), [1, size(N_phase_Mig, 2)]);

freqs_aux = 0:1000;
correctFactTheo = sqrt(1i * freqs_aux/c);

axAbsIndRub = axes(figure);
C = zeros(length(freqEdges), length(absCorrFactEdges));
C(1:end-1, 1:end-1) = Nnorm_abs_Rub;
pcolor(axAbsIndRub, freqEdges, absCorrFactEdges, C');
axAbsIndRub.Title.String = '|\Psi_{ind}| for individual cancellation. (A).';
axAbsIndRub.XLabel.String = 'Frequency (Hz)';
axAbsIndRub.YLabel.String = '|\Psi_{ind}|';
colorbar
axAbsIndRub.NextPlot = 'Add';
plot(axAbsIndRub, freqs_aux, abs(correctFactTheo), 'r', 'LineWidth', 4);

axPhaseIndRub = axes(figure);
C = zeros(length(freqEdges), length(phaseCorrFactEdges));
C(1:end-1, 1:end-1) = Nnorm_phase_Rub;
pcolor(axPhaseIndRub, freqEdges, phaseCorrFactEdges, C');
axPhaseIndRub.Title.String = 'angle(\Psi_{ind}) for individual cancellation. (A).';
axPhaseIndRub.XLabel.String = 'Frequency (Hz)';
axPhaseIndRub.YLabel.String = 'angle(\Psi_{ind}) (º)';
colorbar
axPhaseIndRub.NextPlot = 'Add';
plot(axPhaseIndRub, freqs_aux, rad2deg(angle(correctFactTheo)), 'r', 'LineWidth', 4);

axAbsIndMig = axes(figure);
C = zeros(length(freqEdges), length(absCorrFactEdges));
C(1:end-1, 1:end-1) = Nnorm_abs_Mig;
pcolor(axAbsIndMig, freqEdges, absCorrFactEdges, C');
axAbsIndMig.Title.String = '|\Psi_{ind}| for individual cancellation. (B).';
axAbsIndMig.XLabel.String = 'Frequency (Hz)';
axAbsIndMig.YLabel.String = '|\Psi_{ind}|';
colorbar
axAbsIndMig.NextPlot = 'Add';
plot(axAbsIndMig, freqs_aux, abs(correctFactTheo), 'r', 'LineWidth', 4);

axPhaseIndMig = axes(figure);
C = zeros(length(freqEdges), length(phaseCorrFactEdges));
C(1:end-1, 1:end-1) = Nnorm_phase_Mig;
pcolor(axPhaseIndMig, freqEdges, phaseCorrFactEdges, C');
axPhaseIndMig.Title.String = 'angle(\Psi_{ind}) for individual cancellation. (B).';
axPhaseIndMig.XLabel.String = 'Frequency (Hz)';
axPhaseIndMig.YLabel.String = 'angle(\Psi_{ind}) (º)';
colorbar
axPhaseIndMig.NextPlot = 'Add';
plot(axPhaseIndMig, freqs_aux, rad2deg(angle(correctFactTheo)), 'r', 'LineWidth', 4);

