%% Experiment 9
% Check the difference between Miguel's script and mine
% The conclussion is that the only difference between his processing and
% mine is the generation of the parameter r0 for the attenuation
% calculation, and that he doesn't scale by the loudspeaker separation

%% Preamble
pathSetUp;

imagesPath = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\TFM\Img\';

dataPathName = [globalPath, 'Data\'];
ID = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

%% System set up.

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

% Simulation options
timeDomainActive = true;
fakeTimeProcessing = false;
frequencyDomainActive = false;
automaticLengthModification = false;
predefSignals = true;
saveSignals = true;

%% A) Check that Miguel's scenario parameters produce the same result in both simulations

    %% Simulate with Miguel code
zPos = 1.65;
NSpositions = [0.1 2.3 zPos];
trev = 0; % Tiempo de reverberación

[h_array,WFSpos,tecta,L,mallado_x,mallado_y,h_ad_sources] = SalaGtac(zPos, 1, 1, trev, numSampIR, fs, c, NSpositions);

[filtros_array,an,tn,activo_array]=WFS_DrivingSignals(WFSpos, NSpositions, c, fs);

durSign = 1; % Duration of signal
t = (0:ceil(durSign*fs)-1)/fs;
freqCos = 200;
NSsignal = sin(2*pi*freqCos*t);

[POT_ad,POT_ar]=generamapa(-NSsignal', h_array, filtros_array, activo_array, NSsignal', h_ad_sources(:,1,:,:));

% % Representamos el mapa de presiones
%  dv=50; % Duración de la ventana para promediar al potencia
%  despv=25;  % Desplazamiento del enventanado para realizar el calculo de la
% % potencia
% 
% dibujapot(POT_ad+POT_ar ,L, WFSpos, NSpositions, mallado_x, mallado_y,dv,despv);

    %% Simulate the same scenario but with my code
if ~exist('obj', 'var') || ~isvalid(obj)
    obj = SimulationController;
end
    
roomDim = L;
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

axWFS = axes(figure);
plot(axWFS, t, recWFS_signals_RubScript(1,:), t, d*recWFS_signals_Miguel(1,:), t, recNS_signals_Miguel(1,:))
axWFS.XLabel.String = 'Time (s)';
axWFS.YLabel.String = 'Signal (arbitrary units)';
axWFS.Title.String = 'Received signal from the WFS array';

recNScoef = signal2pulseCoefficientMatrix([0 durSign], freqCos, 1, (recNS_signals)', fs);
recWFScoef = signal2pulseCoefficientMatrix([0 durSign], freqCos, 1, (recWFS_signals_RubScript)', fs);
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

%% B) Simulate the same scenario but with a chirp signal
% The same as previous case, but with a chirp signal. Calculate the
% correction factor

% Chirp signal for the noise source
durSign = 1; % Duration of tone for time processing
t = (0:ceil(durSign*fs)-1)/fs;
NSsignal = chirp(t, 20, durSign, 940);

WFSfilterLength = 2^12;
numSampIR = 2^12;
WFS_AcPath_previously_calculated = true;

SetupParametersScript;
AcousticPathCalculationScript
simulationScript;
rec_time_Chirp = rec_signals;
recNS_time_Chirp = recNS_signals;

% Analysis

% Perform FFT
recNS_freq_Chirp = fft(permute(recNS_time_Chirp, [2 1 3 4]));

    % Rubén attenuation
recWFS_time_Chirp = rec_time_Chirp - recNS_time_Chirp;
recWFS_freq_Chirp = fft(permute(recWFS_time_Chirp, [2 1 3 4])); %/fs;
rec_freq_Chirp = fft(permute(rec_time_Chirp, [2 1 3 4])); %/fs;
corrFact_Chirp = -recNS_freq_Chirp./recWFS_freq_Chirp; % (numSamp x numMicro)

numSamp = size(recNS_time_Chirp, 2);
f = (0:numSamp - 1)*fs/numSamp;
corrFact_phase = rad2deg(angle(corrFact_Chirp));
corrFact_abs = abs(corrFact_Chirp);
ax = axes(figure);
ind = f >=0 & f<1000;
% plot(ax, f(ind), corrFact_abs(ind, :))
plot(ax, f(ind), corrFact_phase(ind, :))

% Conclussion: the phase shift needed is very low for los frequencies, but
% around 40º (45º theoretically) for higher frequencies

%% C) Use multiple noise source positions and compare Miguel's attenuation calculation (parameter r0) and mine

% Positions of the noise source
% Quarter of a circle
numPointsPerArc = 4;
radius = [3.6 4 4.4 4.8];
numArcs = numel(radius);
xOctagon = obj.WFSposition(:, 1);
yOctagon = obj.WFSposition(:, 2);
centreX = (max(xOctagon) + min(xOctagon))/2;
centreY = (max(yOctagon) + min(yOctagon))/2;
y1 = centreY;
y2 = roomDim(2) - centreY;
alphaMax = pi + asin(y1/max(radius)) - deg2rad(1);
alphaMin = pi - asin(y2/max(radius)) + deg2rad(1);
alpha = linspace(alphaMin, alphaMax, numPointsPerArc)';
x = centreX + repmat(radius, numPointsPerQuarter, 1).*repmat(cos(alpha), 1, numArcs);
y = centreY + repmat(radius, numPointsPerQuarter, 1).*repmat(sin(alpha), 1, numArcs);
NSpositions = [x(:), y(:), zPos*ones(numel(x), 1)];

WFS_AcPath_previously_calculated = true;
SetupParametersScript
AcousticPathCalculationScript

attenuationType = 'Miguel';
SetupParametersScript
simulationScript;
rec_signals_Mig = rec_signals;
recNS_signals_Mig = recNS_signals;

attenuationType = 'Ruben';
SetupParametersScript;
simulationScript;
rec_signals_Rub = rec_signals;
recNS_signals_Rub = recNS_signals;

% Analysis
assert(isequal(recNS_signals_Mig, recNS_signals_Rub), 'Signals from NS should be the same')
recNS_time = recNS_signals_Mig;

numSamp = size(recNS_time, 2);
t = (0:numSamp-1)/fs;
f = (0:numSamp-1)*fs/numSamp;
ind = f >= 0 & f <= 1000;
fsel = f(ind);

% Perform FFT
recNS_freq = fft(permute(recNS_time, [2 1 3 4]));
recNS_freq = recNS_freq(ind, :, :, :);

    % Rubén attenuation
recWFS_time_Rub = rec_signals_Rub - recNS_time;
recWFS_freq_Rub = fft(permute(recWFS_time_Rub, [2 1 3 4])); %/fs;
recWFS_freq_Rub = recWFS_freq_Rub(ind, :, :, :);
corrFact_Rub = -recNS_freq./recWFS_freq_Rub;

    % Miguel's attenuation
recWFS_time_Mig = rec_signals_Mig - recNS_time;
recWFS_freq_Mig = fft(permute(recWFS_time_Mig, [2 1 3 4])); %/fs;
recWFS_freq_Mig = recWFS_freq_Mig(ind, :, :, :);
corrFact_Mig = -recNS_freq./recWFS_freq_Mig;

% Graphs
% Histograms
freqEdges = 10:10:1000;
absCorrFactEdges = 0:0.1:4;
phaseCorrFactEdges = -10:90;

freqMat_ind = pointWiseExtend(fsel', corrFact_Rub);
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

% Global cancellation
corrFactGlob_Rub = zeros(length(fsel), numNSpos);
corrFactGlob_Mig = zeros(length(fsel), numNSpos);
for f_ind = 1:length(fsel)
for ns = 1:numNSpos
    nscoef = permute(recNS_freq(f_ind, :, 1, ns), [2, 1]);
    
    wfscoefRub = permute(recWFS_freq_Rub(f_ind, :, 1, ns), [2 1]);
    corrFactGlob_Rub(f_ind, ns) = wfscoefRub\(-nscoef);
    
    wfscoefMig = permute(recWFS_freq_Mig(f_ind, :, 1, ns), [2 1]);
    corrFactGlob_Mig(f_ind, ns) = wfscoefMig\(-nscoef);
end
end

freqMat_ind = pointWiseExtend(fsel', corrFactGlob_Rub);
N_abs_Glob_Rub = histcounts2(freqMat_ind, abs(corrFactGlob_Rub), freqEdges, absCorrFactEdges);
N_phase_Glob_Rub = histcounts2(freqMat_ind, rad2deg(angle(corrFactGlob_Rub)), freqEdges, phaseCorrFactEdges);
N_abs_Glob_Mig = histcounts2(freqMat_ind, abs(corrFactGlob_Mig), freqEdges, absCorrFactEdges);
N_phase_Glob_Mig = histcounts2(freqMat_ind, rad2deg(angle(corrFactGlob_Mig)), freqEdges, phaseCorrFactEdges);

Nnorm_abs_Glob_Rub = N_abs_Glob_Rub./repmat(sum(N_abs_Glob_Rub, 2), [1, size(N_abs_Glob_Rub, 2)]);
Nnorm_phase_Glob_Rub = N_phase_Glob_Rub./repmat(sum(N_phase_Glob_Rub, 2), [1, size(N_phase_Glob_Rub, 2)]);
Nnorm_abs_Glob_Mig = N_abs_Glob_Mig./repmat(sum(N_abs_Glob_Mig, 2), [1, size(N_abs_Glob_Mig, 2)]);
Nnorm_phase_Glob_Mig = N_phase_Glob_Mig./repmat(sum(N_phase_Glob_Mig, 2), [1, size(N_phase_Glob_Mig, 2)]);

axAbsIndRubGlob = axes(figure);
C = zeros(length(freqEdges), length(absCorrFactEdges));
C(1:end-1, 1:end-1) = Nnorm_abs_Glob_Rub;
pcolor(axAbsIndRubGlob, freqEdges, absCorrFactEdges, C');
axAbsIndRubGlob.Title.String = '|\Psi_{ind}| for global cancellation. (A).';
axAbsIndRubGlob.XLabel.String = 'Frequency (Hz)';
axAbsIndRubGlob.YLabel.String = '|\Psi_{ind}|';
colorbar
axAbsIndRubGlob.NextPlot = 'Add';
plot(axAbsIndRubGlob, freqs_aux, abs(correctFactTheo), 'r', 'LineWidth', 4);

axPhaseIndRub = axes(figure);
C = zeros(length(freqEdges), length(phaseCorrFactEdges));
C(1:end-1, 1:end-1) = Nnorm_phase_Glob_Rub;
pcolor(axPhaseIndRub, freqEdges, phaseCorrFactEdges, C');
axPhaseIndRub.Title.String = 'angle(\Psi_{ind}) for global cancellation. (A).';
axPhaseIndRub.XLabel.String = 'Frequency (Hz)';
axPhaseIndRub.YLabel.String = 'angle(\Psi_{ind}) (º)';
colorbar
axPhaseIndRub.NextPlot = 'Add';
plot(axPhaseIndRub, freqs_aux, rad2deg(angle(correctFactTheo)), 'r', 'LineWidth', 4);

axAbsIndMig = axes(figure);
C = zeros(length(freqEdges), length(absCorrFactEdges));
C(1:end-1, 1:end-1) = Nnorm_abs_Glob_Mig;
pcolor(axAbsIndMig, freqEdges, absCorrFactEdges, C');
axAbsIndMig.Title.String = '|\Psi_{ind}| for global cancellation. (B).';
axAbsIndMig.XLabel.String = 'Frequency (Hz)';
axAbsIndMig.YLabel.String = '|\Psi_{ind}|';
colorbar
axAbsIndMig.NextPlot = 'Add';
plot(axAbsIndMig, freqs_aux, abs(correctFactTheo), 'r', 'LineWidth', 4);

axPhaseIndMig = axes(figure);
C = zeros(length(freqEdges), length(phaseCorrFactEdges));
C(1:end-1, 1:end-1) = Nnorm_phase_Glob_Mig;
pcolor(axPhaseIndMig, freqEdges, phaseCorrFactEdges, C');
axPhaseIndMig.Title.String = 'angle(\Psi_{ind}) for global cancellation. (B).';
axPhaseIndMig.XLabel.String = 'Frequency (Hz)';
axPhaseIndMig.YLabel.String = 'angle(\Psi_{ind}) (º)';
colorbar
axPhaseIndMig.NextPlot = 'Add';
plot(axPhaseIndMig, freqs_aux, rad2deg(angle(correctFactTheo)), 'r', 'LineWidth', 4);

%% D) Multiple NS positions not limited by room dimension, time and frequency processing
% Positions of the noise source
% Quarter of a circle
numPointsPerArc = 6;
radius = [5 7.5 10];
numArcs = numel(radius);
xOctagon = obj.WFSposition(:, 1);
yOctagon = obj.WFSposition(:, 2);
centreX = (max(xOctagon) + min(xOctagon))/2;
centreY = (max(yOctagon) + min(yOctagon))/2;
alphaMax = pi/2;
alphaMin = 0;
alpha = linspace(alphaMin, alphaMax, numPointsPerArc)';
x = centreX + repmat(radius, numPointsPerQuarter, 1).*repmat(cos(alpha), 1, numArcs);
y = centreY + repmat(radius, numPointsPerQuarter, 1).*repmat(sin(alpha), 1, numArcs);
NSpositions = [x(:), y(:), zPos*ones(numel(x), 1)];

WFS_AcPath_previously_calculated = true;
SetupParametersScript
AcousticPathCalculationScript

attenuationType = 'Miguel';
SetupParametersScript
simulationScript;
rec_signals_Mig = rec_signals;
recNS_signals_Mig = recNS_signals;

attenuationType = 'Ruben';
SetupParametersScript;
simulationScript;
rec_signals_Rub = rec_signals;
recNS_signals_Rub = recNS_signals;


%% E) Use even Miguel's script parameters and processing, but with a chirp signal
% Be aware that fs will change to the value set in the script

%%%%%%%%%%%%%%%%%%%% ejemplocarwfs.m code
fs = 44100;
c=340;
maxDist = 20; % meteres
numSampIR = ceil(maxDist/c*fs);

fte = [0.1 2.3 1.65;0.1 1.3 1.65;0.1 3.3 1.65];
%fte=[alt(1,93)-0.001 alt(2,93) 1.65;2,2.3,1.65];
%fte=alt(:,93);
fuente_ruido=1;

% Definimos la sala y las fuentes de ruido
disp('Calculando las respuestas del sistema acústico a modelar...')
% [h_array,alt,tecta,L,mallado_x,mallado_y,h_ad_sources] = SalaGtac(1.65,0.05,0.05,0,250,fs,c,fte,1);
% save salasinrev h_array h_ad_sources alt tecta fs fte L mallado_x mallado_y;
% 
% 
% [h_array,alt,tecta,L,mallado_x,mallado_y,h_ad_sources] = SalaGtac(1.65,0.2,0.2,0.15,250,fs,c,fte);
% save salaconrev2 h_array h_ad_sources alt tecta fs fte L mallado_x mallado_y;
[h_array,alt,tecta,L,mallado_x,mallado_y,h_ad_sources] = SalaGtac(1.65,1, 1,0,numSampIR,fs,c,fte(fuente_ruido,:));
%save salasinr3 h_array h_ad_sources alt tecta fs fte L mallado_x mallado_y;
% Calculamos las funciones directoras del array definido por 'alt' para sintentizar una fuente en la posición fte. 

disp('Calculando la configuración del array...')
[filtros_array,an,tn,activo_array]=WFS_DrivingSignals(alt,fte(fuente_ruido,:),c,fs);

% Definimoa la señal de ruido
durSign = 1; % Duration of tone for time processing
t = (0:ceil(durSign*fs)-1)/fs;
NSsignal = chirp(t, 20, durSign, 940)';

[POT_ad,POT_ar]=generamapa(-NSsignal,h_array,filtros_array,activo_array,NSsignal,h_ad_sources(:,fuente_ruido,:,:));
%%%%%%%%%%%%%%%%%%%%

% For comparative analysis
eval(['POT_ad_', num2str(fs), ' = POT_ad;'])
eval(['POT_ar_', num2str(fs), ' = POT_ar;'])
eval(['h_ad_sources_', num2str(fs), ' = h_ad_sources;'])
eval(['h_array_', num2str(fs), ' = h_array;'])
eval(['NSsignal_', num2str(fs), ' = NSsignal;'])

% Analysis
% Perform FFT
recNS_freq = fft(POT_ad);

    % Rubén attenuation
recWFS_freq = fft(POT_ar);
corrFact = -recNS_freq./recWFS_freq; % (numSamp x numMicro)

numSamp = size(POT_ad, 1);
f = (0:numSamp - 1)*fs/numSamp;
ind = f >=0 & f<1000;
corrFactRed = corrFact(ind, 1:3:size(corrFact, 2));
corrFact_phase = rad2deg(angle(corrFactRed));
corrFact_abs = abs(corrFactRed);
ax = axes(figure);
plot(ax, f(ind), corrFact_abs(ind, 1))
plot(ax, f(ind), corrFact_phase(ind, :))

% Comparative analysis: different fs
numSamp8000 = size(POT_ad_8000, 1);
numSamp44100 = size(POT_ad_44100, 1);
t8000 = (0:numSamp8000 - 1)/8000;
t44100 = (0:numSamp44100 - 1)/44100;
f8000 = (0:numSamp8000 - 1)*8000/numSamp8000;
f44100 = (0:numSamp44100 - 1)*44100/numSamp44100;
indMicro = 10;

%
recNS_freq_8000 = fft(POT_ad_8000);
recNS_freq_44100 = fft(POT_ad_44100);

recNS_freq_8000_phase = rad2deg(unwrap(angle(recNS_freq_8000(:, indMicro))));
recNS_freq_44100_phase = rad2deg(unwrap(angle(recNS_freq_44100(:, indMicro))));
recNS_freq_44100_phase_adapt = interp1(f44100, recNS_freq_44100_phase, f8000');

axNS = axes(figure);
% plot(axNS, f8000, abs(recNS_freq_8000(:, indMicro)*(44100/8000)), f44100, abs(recNS_freq_44100(:, indMicro)))
plot(axNS, f8000, recNS_freq_8000_phase, f44100, recNS_freq_44100_phase)
axNS.XLim = [0, 1000];

axNSdif = axes(figure);
plot(axNSdif, f8000, recNS_freq_8000_phase - recNS_freq_44100_phase_adapt);
axNSdif.XLim = [0, 1000];

axNStime = axes(figure);
plot(axNStime, t8000, POT_ad_8000(:, indMicro), t44100, POT_ad_44100(:, indMicro));

[xInd, yInd] = ind2sub([length(mallado_x), length(mallado_y)], indMicro);
tAcPath_8000 = (0:size(h_ad_sources_8000, 1)-1)/8000;
tAcPath_44100 = (0:size(h_ad_sources_44100, 1)-1)/44100;
axAcPath = axes(figure);
plot(axAcPath, tAcPath_8000, h_ad_sources_8000(:, 1, xInd, yInd),...
    tAcPath_44100, h_ad_sources_44100(:, 1, xInd, yInd))
% It seems the impulse response is actually good
% Filter the noise source signal to check if it produces some delay
filt8000 = filter(h_ad_sources_8000(:, 1, xInd, yInd), 1, NSsignal_8000);
filt44100 = filter(h_ad_sources_44100(:, 1, xInd, yInd), 1, NSsignal_44100);
axFilt = axes(figure);
plot(axFilt, t8000, filt8000, t44100, filt44100)

%

%
recWFS_freq_8000 = fft(POT_ar_8000);
recWFS_freq_44100 = fft(POT_ar_44100);

recWFS_freq_8000_phase = rad2deg(unwrap(angle(recWFS_freq_8000(:, indMicro))));
recWFS_freq_44100_phase = rad2deg(unwrap(angle(recWFS_freq_44100(:, indMicro))));

recWFS_freq_44100_phase_adapt = interp1(f44100, recWFS_freq_44100_phase, f8000');

axWFS = axes(figure);
% plot(axWFS, f8000, abs(recWFS_freq_8000(:, indMicro)*(44100/8000)), f44100, abs(recWFS_freq_44100(:, indMicro)))
plot(axWFS, f8000, recWFS_freq_8000_phase, f44100, recWFS_freq_44100_phase)
axWFS.XLim = [0, 1000];

axWFSSdif = axes(figure);
plot(axWFSSdif, f8000, recWFS_freq_8000_phase - recWFS_freq_44100_phase_adapt);
axWFSSdif.XLim = [0, 1000];

axWFStime = axes(figure);
plot(axWFStime, t8000, POT_ar_8000(:, indMicro), t44100, POT_ar_44100(:, indMicro))
% There is a delay between signals that produces a linear phase shift with
% frequency
%

%
corrFact_8000 = -recNS_freq_8000./recWFS_freq_8000;
corrFact_44100 = -recNS_freq_44100./recWFS_freq_44100;

corrFact_8000_phase = rad2deg(unwrap(angle(corrFact_8000(:, indMicro))));
corrFact_44100_phase = rad2deg(unwrap(angle(corrFact_44100(:, indMicro))));

corrFact_44100_phase_adapt = interp1(f44100, corrFact_44100_phase, f8000');

axCorr = axes(figure);
% plot(axCorr, f8000, abs(corrFact_8000(:, indMicro)), f44100, abs(corrFact_44100(:, indMicro)))
% axCorr.XLim = [0, 1000];
plot(axCorr, f8000, corrFact_8000_phase, f44100, corrFact_44100_phase)
axCorr.XLim = [0, 1000];

axCorrDiff = axes(figure);
plot(axCorrDiff, f8000, corrFact_8000_phase - corrFact_44100_phase_adapt)
axCorrDiff.XLim = [0, 1000];
%

% Conclussion: different fs (sampling frequency) result in very different
% phase shifts. This happens because for difference sampling frequencies,
% the delay of the signals (at least the received WFS ones) is different, hence there is a linear phase
% shift. Why is the delay different?
% I have a hypothesis.
% Since sampling introduces a quantification error in the time domain, the
% mangitude of that error depends on the sampling frequency. For example,
% for a sampling frequency of 10Hz, and a delay of 1.9 seconds, one must
% choose to put the delta at 0.1 or 0.2, depending if one uses the function
% floor, ceil or round to discretize the positions. So, the maximum
% precission error is 1/fs. At a given frequency f, the phase shift that this
% delay introduces is (1/fs)*f*360 degrees.
% So, if we want a maximum error of E degrees:
% (1/fs)*f*360 < E --> fs > f*360/E.
% For a maximum error of 10º, and a maximum frequency of 940Hz (aliasing frequency),
% fs > 33840 Hz. So, a sampling frequency of 44100 is not only adequate,
% but it can even be not hight enough for good precissions.