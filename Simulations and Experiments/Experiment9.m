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
plot(axNS, t, recNS_signals(15,:), t, recNS_signals_Miguel(15,:))
axNS.XLabel.String = 'Time (s)';
axNS.YLabel.String = 'Signal (arbitrary units)';
axNS.Title.String = 'Received signal from the noise source';

dif = recNS_signals(:, 1:size(POT_ad, 1)) - recNS_signals_Miguel(:, 1:size(POT_ad, 1));
max(abs(dif), [], 2)

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
obj.domain = 'time';
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
x = centreX + repmat(radius, numPointsPerArc, 1).*repmat(cos(alpha), 1, numArcs);
y = centreY + repmat(radius, numPointsPerArc, 1).*repmat(sin(alpha), 1, numArcs);
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

    % Rubén attenuation
recWFS_time_Rub = rec_signals_Rub - recNS_time;
[corrFact_Rub, corrFactGlob_Rub, fsel, axAbsIndRub, axPhaseIndRub, axAbsGlobRub, axPhaseGlobRub] = calculateTransferFunction(recWFS_time_Rub, -recNS_time, fs, {'minimizDim', 'time', 'reverb', 'NSpos'});

corrFactTheo = sqrt(1i*fsel/c);

axAbsIndRub.NextPlot = 'Add';
plot(axAbsIndRub, fsel, abs(corrFactTheo), 'r', 'LineWidth', 3);
axAbsGlobRub.NextPlot = 'Add';
plot(axAbsGlobRub, fsel, abs(corrFactTheo), 'r', 'LineWidth', 3);

% printfig(axAbsIndRub.Parent, imagesPath, 'Experiment9_CorrFactAbsIndChirpRubAtten', 'eps');
% printfig(axPhaseIndRub.Parent, imagesPath, 'Experiment9_CorrFactAbsIndChirpRubAtten', 'eps');

    % Miguel's attenuation
recWFS_time_Mig = rec_signals_Mig - recNS_time;
[~, ~, ~, axAbsIndMig, axPhaseIndMig, axAbsGlobMig, axPhaseGlobMib] = calculateTransferFunction(recWFS_time_Mig, -recNS_time, fs, {'minimizDim', 'time', 'reverb', 'NSpos'});

axPhaseIndMig.YLim = [0 50];
axPhaseGlobMib.YLim = [0 50];
axAbsIndMig.YLim = [0 3];
axAbsGlobMig.YLim = [0 3];

axAbsIndMig.NextPlot = 'Add';
plot(axAbsIndMig, fsel, abs(corrFactTheo), 'r', 'LineWidth', 3);
axAbsGlobMig.NextPlot = 'Add';
plot(axAbsGlobMig, fsel, abs(corrFactTheo), 'r', 'LineWidth', 3);

% printfig(axAbsIndMig.Parent, imagesPath, 'Experiment9_CorrFactAbsIndChirpMigAtten', 'eps');
% printfig(axPhaseIndMig.Parent, imagesPath, 'Experiment9_CorrFactPhaseIndChirpMigAtten', 'eps');

% % Graphs
% domain = [0 1];
% visualObj = animation({fsel, 1:numMicro, 1:numReverbTime, 1:numNSpos},...
%     {abs(corrFact_Rub)}, {'Frequency', 'Microphone', 'Reverb. Time index', 'NS position index'}, {'Cancellation'}, [], []);

% Now simulate in the frequency domain
timeDomainActive = false;
fakeTimeProcessing = false;
frequencyDomainActive = true;
saveSignals = false;
freqs = fsel(1:4:end);

SetupParametersScript
AcousticPathCalculationScript

attenuationType = 'Miguel';
SetupParametersScript
simulationScript;
sMig = s;

attenuationType = 'Ruben';
SetupParametersScript;
simulationScript;
sRub = s;

% Calculate individual and global correction factors
sMigFreq = sMig(2, :, :, :);
[sMigFreq, corrFactIndMig, corrFactGlobalMig] = SimulationController.addCancellationParametersToStructure(sMigFreq);

sRubFreq = sRub(2, :, :, :);
[sRubFreq, corrFactIndRub, corrFactGlobalRub] = SimulationController.addCancellationParametersToStructure(sRubFreq);

% Hisgoram
absCorrFactEdges = 0:0.1:4;
phaseCorrFactEdges = -10:90;

% Miguel
axAbsIndMigFreq = histogram2D( abs(corrFactIndMig), 2, freqs, [], absCorrFactEdges );
axPhaseIndMigFreq = histogram2D( rad2deg(angle(corrFactIndMig)), 2, freqs, [], phaseCorrFactEdges );
axPhaseIndMigFreq.YLim = [0 50];

axAbsGlobMigFreq = histogram2D( abs(corrFactGlobalMig), 2, freqs, [], absCorrFactEdges );
axPhaseGlobMigFreq = histogram2D( rad2deg(angle(corrFactGlobalMig)), 2, freqs, [], phaseCorrFactEdges );
axPhaseGlobMigFreq.YLim = [0 50];

axAbsIndMigFreq.Parent.Name = 'CorrFact Ind Mig Abs. Freq processing.';
axPhaseIndMigFreq.Parent.Name = 'CorrFact Ind Mig Phase. Freq processing.';
axAbsGlobMigFreq.Parent.Name = 'CorrFact Glob Mig Abs. Freq processing.';
axPhaseGlobMigFreq.Parent.Name = 'CorrFact Glob Mig Phase. Freq processing.';

axAbsIndMigFreq.NextPlot = 'Add';
plot(axAbsIndMigFreq, fsel, abs(corrFactTheo), 'r', 'LineWidth', 3);
axAbsGlobMigFreq.NextPlot = 'Add';
plot(axAbsGlobMigFreq, fsel, abs(corrFactTheo), 'r', 'LineWidth', 3);

% printfig(axAbsIndMigFreq.Parent, imagesPath, 'Experiment9_CorrFactAbsIndFreqMigAtten', 'eps');
% printfig(axPhaseIndMigFreq.Parent, imagesPath, 'Experiment9_CorrFactAbsIndFreqMigAtten', 'eps');

% Rubén
axAbsIndRubFreq = histogram2D( abs(corrFactIndRub), 2, freqs, [], absCorrFactEdges );
axPhaseIndRubFreq = histogram2D( rad2deg(angle(corrFactIndRub)), 2, freqs, [], phaseCorrFactEdges );
axPhaseIndRubFreq.YLim = [0 50];

axAbsGlobRubFreq = histogram2D( abs(corrFactGlobalRub), 2, freqs, [], absCorrFactEdges );
axPhaseGlobRubFreq = histogram2D( rad2deg(angle(corrFactGlobalRub)), 2, freqs, [], phaseCorrFactEdges );
axPhaseGlobRubFreq.YLim = [0 50];

axAbsIndRubFreq.Parent.Name = 'CorrFact Ind Rub Abs. Freq processing.';
axPhaseIndRubFreq.Parent.Name = 'CorrFact Ind Rub Phase. Freq processing.';
axAbsGlobRubFreq.Parent.Name = 'CorrFact Glob Rub Abs. Freq processing.';
axPhaseGlobRubFreq.Parent.Name = 'CorrFact Glob Rub Phase. Freq processing.';

axAbsGlobRubFreq.NextPlot = 'Add';
plot(axAbsGlobRubFreq, fsel, abs(corrFactTheo), 'r', 'LineWidth', 3);
axAbsIndRubFreq.NextPlot = 'Add';
plot(axAbsIndRubFreq, fsel, abs(corrFactTheo), 'r', 'LineWidth', 3);

% SVG scenario
viewBox = [-WFSarrayOffset(1) -WFSarrayOffset(2) roomDim(1) roomDim(2)];
centreX = (max(obj.WFSposition(:, 1)) + min(obj.WFSposition(:, 1)))/2;
centreY = (max(obj.WFSposition(:, 2)) + min(obj.WFSposition(:, 2)))/2;
NSangles = atan2d(centreY - NSpositions(:,2), centreX - NSpositions(:,1));

objSVG = SVGdrawer('viewBox', viewBox, 'NSpositions', NSpositions,...
    'NSangles', NSangles, 'microSymbol', 'dot', 'microSize', 0.05,...
    'microPositions', recPositions, 'WFSpositions', obj.WFSposition(:, [1 2]), 'WFSangles', tecta - 90);

name = 'Experiment9_NSinsideChamberScheme';
% objSVG.drawSVG([imagesPath, name, '.svg']);
% 
% currentFolder = pwd;
% cd(imagesPath); % Needed for inkscape to link svg files properly
% system(['inkscape -z "', imagesPath, name, '.svg" --export-pdf="', imagesPath, name, '.pdf"'])
% cd(currentFolder)

%% D) Multiple NS positions not limited by room dimension. Use time and frequency processing
% Positions of the noise source
% Arc of circle
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
x = centreX + repmat(radius, numPointsPerArc, 1).*repmat(cos(alpha), 1, numArcs);
y = centreY + repmat(radius, numPointsPerArc, 1).*repmat(sin(alpha), 1, numArcs);
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

    % Rubén attenuation
recWFS_time_Rub = rec_signals_Rub - recNS_time;
[corrFact_Rub, corrFactGlob_Rub, fsel, axAbsIndRub, axPhaseIndRub, axAbsGlobRub, axPhaseGlobRub] = calculateTransferFunction(recWFS_time_Rub, -recNS_time, fs, {'minimizDim', 'time', 'reverb', 'NSpos'});
axAbsIndRub.Parent.Name = 'CorrFact Ind Rub Abs. Chirp.';
axPhaseIndRub.Parent.Name = 'CorrFact Ind Rub Phase. Chirp';
axAbsGlobRub.Parent.Name = 'CorrFact Glob Rub Abs. Chirp';
axPhaseGlobRub.Parent.Name = 'CorrFact Glob Rub Phase. Chirp';
axPhaseIndRub.YLim = [0 50];
axPhaseGlobRub.YLim = [0 50];

    % Miguel's attenuation
recWFS_time_Mig = rec_signals_Mig - recNS_time;
[corrFact_Mig, corrFactGlob_Mig, ~, axAbsIndMig, axPhaseIndMig, axAbsGlobMig, axPhaseGlobMib] = calculateTransferFunction(recWFS_time_Mig, -recNS_time, fs, {'minimizDim', 'time', 'reverb', 'NSpos'});
axAbsIndMig.Parent.Name = 'CorrFact Ind Mig Abs. Chirp';
axPhaseIndMig.Parent.Name = 'CorrFact Ind Mig Phase. Chirp';
axAbsGlobMig.Parent.Name = 'CorrFact Glob Mig Abs. Chirp';
axPhaseGlobMig.Parent.Name = 'CorrFact Glob Mig Phase. Chirp';
axPhaseIndMig.YLim = [0 50];
axPhaseGlobMig.YLim = [0 50];

% Now simulate in the frequency domain
timeDomainActive = false;
fakeTimeProcessing = false;
frequencyDomainActive = true;
saveSignals = false;
freqs = fsel;

SetupParametersScript
AcousticPathCalculationScript

attenuationType = 'Miguel';
SetupParametersScript
simulationScript;
sMig = s;

% dimensionOrder = {'domain', 'frequency', 'NSposition', 'ReverbTime'};
% dataStruct = struct();
% dataStruct.domain = {'frequency', 'time'};
% dataStruct.frequency = freqs;
% dataStruct.NSposition = NSpositions;
% dataStruct.ReverbTime = beta;
% dataStruct.dimensionOrder = dimensionOrder;
% dataStruct.data = sMig;
% dataStruct.recPositions = recPositions;
% 
% save([dataPathName, 'Experiment9sMig', ID, '.mat'], 'dataStruct')


attenuationType = 'Ruben';
SetupParametersScript;
simulationScript;
sRub = s;

% dimensionOrder = {'domain', 'frequency', 'NSposition', 'ReverbTime'};
% dataStruct = struct();
% dataStruct.domain = {'frequency', 'time'};
% dataStruct.frequency = freqs;
% dataStruct.NSposition = NSpositions;
% dataStruct.ReverbTime = beta;
% dataStruct.dimensionOrder = dimensionOrder;
% dataStruct.data = sMig;
% dataStruct.recPositions = recPositions;
% 
% save([dataPathName, 'Experiment9sRub', ID, '.mat'], 'dataStruct')

% Calculate individual and global correction factors
sMigFreq = sMig(2, :, :, :);
[sMigFreq, corrFactIndMig, corrFactGlobalMig] = SimulationController.addCancellationParametersToStructure(sMigFreq);

sRubFreq = sRub(2, :, :, :);
[sRubFreq, corrFactIndRub, corrFactGlobalRub] = SimulationController.addCancellationParametersToStructure(sRubFreq);

% Hisgoram
absCorrFactEdges = 0:0.1:4;
phaseCorrFactEdges = -10:90;

axAbsIndMig = histogram2D( abs(corrFactIndMig), 2, freqs, [], absCorrFactEdges );
axPhaseIndMig = histogram2D( rad2deg(angle(corrFactIndMig)), 2, freqs, [], phaseCorrFactEdges );
axPhaseIndMig.YLim = [0 50];

axAbsGlobMig = histogram2D( abs(corrFactGlobalMig), 2, freqs, [], absCorrFactEdges );
axPhaseGlobMig = histogram2D( rad2deg(angle(corrFactGlobalMig)), 2, freqs, [], phaseCorrFactEdges );
axPhaseGlobMig.YLim = [0 50];

axAbsIndMig.Parent.Name = 'CorrFact Ind Mig Abs. Freq processing.';
axPhaseIndMig.Parent.Name = 'CorrFact Ind Mig Phase. Freq processing.';
axAbsGlobMig.Parent.Name = 'CorrFact Glob Mig Abs. Freq processing.';
axPhaseGlobMig.Parent.Name = 'CorrFact Glob Mig Phase. Freq processing.';

axAbsIndRub = histogram2D( abs(corrFactIndRub), 2, freqs, [], absCorrFactEdges );
axPhaseIndRub = histogram2D( rad2deg(angle(corrFactIndRub)), 2, freqs, [], phaseCorrFactEdges );
axPhaseInd.YLim = [0 50];

axAbsGlobRub = histogram2D( abs(corrFactGlobalRub), 2, freqs, [], absCorrFactEdges );
axPhaseGlobRub = histogram2D( rad2deg(angle(corrFactGlobalRub)), 2, freqs, [], phaseCorrFactEdges );
axPhaseGlob.YLim = [0 50];

axAbsIndRub.Parent.Name = 'CorrFact Ind Rub Abs. Freq processing.';
axPhaseIndRub.Parent.Name = 'CorrFact Ind Rub Phase. Freq processing.';
axAbsGlobRub.Parent.Name = 'CorrFact Glob Rub Abs. Freq processing.';
axPhaseGlobRub.Parent.Name = 'CorrFact Glob Rub Phase. Freq processing.';


%% E) Use even Miguel's script parameters and processing, but with a chirp signal
% Be aware that fs will change to the value set in the script

%%%%%%%%%%%%%%%%%%%% ejemplocarwfs.m code
fs = 44100;
c=340;
maxDist = 20; % meteres
numSampIR = ceil(maxDist/c*fs);

fte = [0.1 2.3 1.65;0.1 1.3 1.65;0.1 3.3 1.65];

% Definimoa la señal de ruido
durSign = 1; % Duration of tone for time processing
t = (0:ceil(durSign*fs)-1)/fs;
NSsignal = chirp(t, 20, durSign, 940);

disp('Calculando las respuestas del sistema acústico a modelar...')
[h_array,alt,tecta,L,mallado_x,mallado_y,h_ad_sources] = SalaGtac(1.65,1, 1,0,numSampIR,fs,c,fte);

numNS = size(fte, 1);
numMicro = size(h_array, 3) * size(h_array, 4);

recNSsignalsMig = zeros(length(t), numMicro, numNS);
recWFSsignalsMig = zeros(length(t), numMicro, numNS);
for fuente_ruido = 1:3
% Definimos la sala y las fuentes de ruido
disp('Calculando la configuración del array...')
[filtros_array,an,tn,activo_array]=WFS_DrivingSignals(alt,fte(fuente_ruido,:),c,fs);

% Calculamos señales recibidas
disp('Calculando las señales en los puntos de control...')
[POT_ad,POT_ar]=generamapa(-NSsignal',h_array,filtros_array,activo_array,NSsignal',h_ad_sources(:,fuente_ruido,:,:));

recNSsignalsMig(:, :, fuente_ruido) = POT_ad;
recWFSsignalsMig(:, :, fuente_ruido) = POT_ar;
end
%%%%%%%%%%%%%%%%%%%%

% Check the result is the same as with my code
% To be completed

% Analysis
d = 0.18;
[corrFact, corrFactGlob, fsel, axAbsInd, axPhaseInd, axAbsGlob, axPhaseGlob] =...
    calculateTransferFunction(recWFSsignalsMig*d, -recNSsignalsMig, fs, {'time', 'minimizDim', 'nsPos'});

corrFactTheo = sqrt(1i*fsel/c);

axAbsInd.YLim = [0, 3];
axAbsInd.NextPlot = 'Add';
plot(axAbsInd, fsel, abs(corrFactTheo), 'r', 'LineWidth', 3);

axPhaseInd.YLim = [0, 50];

axAbsGlob.YLim = [0, 3];
axAbsGlob.NextPlot = 'Add';
plot(axAbsGlob, fsel, abs(corrFactTheo), 'r', 'LineWidth', 3);

axPhaseGlob.YLim = [0, 50];

% printfig(axAbsInd.Parent, imagesPath, 'Experiment9_MigScriptCorrFactAbsInd', 'eps');
% printfig(axPhaseInd.Parent, imagesPath, 'Experiment9_MigScriptCorrFactPhaseInd', 'eps');

% % Other analysis
% corrFactResh = reshape(corrFact, [length(fsel), numMicro*numNS]);
% ax = axes(figure);
% lineVec = plot(ax, fsel, abs(corrFactResh));
% 
% % Each noise source with a line color
% cmap = lines;
% for ns = 1:numNS
%     inds = sub2ind([numMicro, numNS], 1:numMicro, ns*ones(1, numMicro));
%     colors = repmat({cmap(ns, :)}, [numMicro, 1]);
%     [lineVec(inds).Color] = colors{:};
% end
% % Nothing noticeable
% 
% % Each control point (microphone position) with a line color
% cmap = lines;
% for micro = 1:numMicro
%     inds = sub2ind([numMicro, numNS], micro*ones(1, numNS), 1:numNS);
%     colors = repmat({cmap(micro, :)}, [numNS, 1]);
%     [lineVec(inds).Color] = colors{:};
% end
% % The lines for the same microphone position tend to get clutterend, maybe,
% % I'm not sure

% Use my code to see if the result is the same
% Important!! First, execute the third section of this m-file. It sets a lot of
% parameters

% Then
if ~exist('obj', 'var') || ~isvalid(obj)
    obj = SimulationController;
end
    
zPos = 1.65;
roomDim = L;
predefRoomDim = true;
[X, Y] = ndgrid(mallado_x, mallado_y);
recPositions = [X(:), Y(:), zPos*ones(numel(X), 1)];
NSpositions = fte;

WFSposition = [alt', zPos*ones(96, 1)];
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

recNSsignalsRub = permute(recNS_signals(:, 1:length(t), 1, :), [2 1 4 3]);
ax = axes(figure);
plot(ax, t, recNSsignalsRub(:, 15, 2), t, recNSsignalsMig(:, 15, 2))

dif = recNSsignalsRub - recNSsignalsMig;
max(abs(dif(:)))

recSignalsRub = permute(rec_signals(:, 1:length(t), 1, :), [2 1 4 3]);
recWFSsignalsRub = recSignalsRub - recNSsignalsRub;
ax = axes(figure);
plot(ax, t, recWFSsignalsRub(:, 1, 3), t, d*recWFSsignalsMig(:, 1, 3))

dif = recWFSsignalsRub - d*recWFSsignalsMig;
max(abs(dif(:)))

% They are different.
% This is due to de fact that the WFS filters are different sometimes by
% one sample. This is because I use this line to calculate the position of
% the delta: 
% indDelta = floor(delays*obj.Fs) + 1;
% And he uses (traduced to my variable names and my algorithm):
% indDelta = round(delays*obj.Fs);
% Change it and the result will be the same.

% % SVG scenario
% viewBox = [0 0 L(1) L(2)];
% centreX = (max(alt(1, :)) + min(alt(1, :)))/2;
% centreY = (max(alt(2, :)) + min(alt(2, :)))/2;
% NSangles = atan2d(centreY - fte(:,2), centreX - fte(:,1));
% 
% [X, Y] = ndgrid(mallado_x, mallado_y);
% recPositions = [X(:), Y(:)];
% 
% objSVG = SVGdrawer('viewBox', viewBox, 'NSpositions', fte,...
%     'NSangles', NSangles, 'microSymbol', 'dot', 'microSize', 0.05,...
%     'microPositions', recPositions, 'WFSpositions', alt', 'WFSangles', tecta - 90);
% 
% name = 'Experiment9_MigScriptScheme';
% % objSVG.drawSVG([imagesPath, name, '.svg']);
% % 
% % currentFolder = pwd;
% % cd(imagesPath); % Needed for inkscape to link svg files properly
% % system(['inkscape -z "', imagesPath, name, '.svg" --export-pdf="', imagesPath, name, '.pdf"'])
% % cd(currentFolder)

%% F) Comparative analysis between sampling frequencies

%%%%%%%%%%%%%%%%%%%% ejemplocarwfs.m code
fs = 44100;
c=340;
maxDist = 20; % meteres
numSampIR = ceil(maxDist/c*fs)1;

fte = [0.1 2.3 1.65;0.1 1.3 1.65;0.1 3.3 1.65];

% Definimoa la señal de ruido
durSign = 1; % Duration of tone for time processing
t = (0:ceil(durSign*fs)-1)/fs;
NSsignal = chirp(t, 20, durSign, 940);

disp('Calculando las respuestas del sistema acústico a modelar...')
[h_array,alt,tecta,L,mallado_x,mallado_y,h_ad_sources] = SalaGtac(1.65,1, 1,0,numSampIR,fs,c,fte);

fuente_ruido = 1;
% Definimos la sala y las fuentes de ruido
disp('Calculando la configuración del array...')
[filtros_array,an,tn,activo_array]=WFS_DrivingSignals(alt,fte(fuente_ruido,:),c,fs);

% Calculamos señales recibidas
disp('Calculando las señales en los puntos de control...')
[POT_ad,POT_ar]=generamapa(-NSsignal',h_array,filtros_array,activo_array,NSsignal',h_ad_sources(:,fuente_ruido,:,:));
%%%%%%%%%%%%%%%%%%%%

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

%% G) Apply theoretical correction factor optimized for one centered frequency
% Noise source positions inside chamber, so the code is very similar to C)
% We don't compare Miguel and Rubén attenuation, it is not relevant to this
% experiment. Just choose one.
% Execute the System set up section first. Then perform a simulation with
% ejemplocarwfs to generate the scenario variables that will be emulated

fopt = 500; % Frequency at which we are going to optimize (Hz)
corrFactTheoOpt = sqrt(1i*fopt/c);
delay = (2*pi - angle(corrFactTheoOpt))/(2*pi*fopt);
ind = floor(delay*fs) + 1;
corrFilter = zeros(ind, 1);
corrFilter(ind) = abs(corrFactTheoOpt);

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

% Chirp signal for the noise source
durSign = 1; % Duration of tone for time processing
t = (0:ceil(durSign*fs)-1)/fs;
NSsignal = chirp(t, 20, durSign, 940);

maxDist = 20;
numSampIR = 2^nextpow2(floor(maxDist/c*fs) + 1);

WFSarrayOffset = [0, 0, 0]; % [x, y, z] coordinates. Useful for generating acoustic path IR.

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
x = centreX + repmat(radius, numPointsPerArc, 1).*repmat(cos(alpha), 1, numArcs);
y = centreY + repmat(radius, numPointsPerArc, 1).*repmat(sin(alpha), 1, numArcs);
NSpositions = [x(:), y(:), zPos*ones(numel(x), 1)];

WFS_AcPath_previously_calculated = true;
attenuationType = 'Miguel';

SetupParametersScript
AcousticPathCalculationScript

simulationScript;

% Analysis
recWFS_signals = rec_signals - recNS_signals;
recWFS_signals_Corr = filter(corrFilter, 1, recWFS_signals, [], 2);

numSamp = size(rec_signals, 2);
t = (0:numSamp - 1)/fs;
ax = axes(figure);
plot(ax, t, recNS_signals(1,:,1,1), t, recWFS_signals_Corr(1,:,1,1) + recNS_signals(1, :, 1, 1));
legend(ax, 'NS', 'Total')
ax.XLabel.String = 'Time (s)';

% printfig(ax.Parent, imagesPath, 'Experiment9_cancOpt460HzTime', 'eps')

% Calculate cancellation
freqsFFT = (0:numSamp - 1)*fs/numSamp;
sel = freqsFFT >= 0 & freqsFFT <= 1000;
freqsFFTsel = freqsFFT(sel);

recNS_freq = fft(recNS_signals, [], 2);
rec_freq = fft(rec_signals, [], 2);
recCorr_freq = fft(recNS_signals + recWFS_signals_Corr, [], 2);

recNS_freq = recNS_freq(:, sel, :, :);
rec_freq = rec_freq(:, sel, :, :);
recCorr_freq = recCorr_freq(:, sel, :, :);

canc = 10*log10(abs(rec_freq./recNS_freq));
cancCorr = 10*log10(abs(recCorr_freq./recNS_freq));

axCanc = histogram2D( canc, 2, freqsFFTsel, [], [] );
axCancCorr = histogram2D( cancCorr, 2, freqsFFTsel, [], [] );
colorbar(axCancCorr)
axCancCorr.YLabel.String = 'Cancellation (dB)';
axCancCorr.XLabel.String = 'Frequency (º)';

% printfig(axCancCorr.Parent, imagesPath, 'Experiment9_cancOpt460Hz', 'eps')

%% H) Diferentes coeficientes de reflexión con corrección frecuencial
% Filter variables for the time WFS filter.
% Constants
c = 340; % Sound velocity (m/s)
fs = 44100; % Sample frequency (samples/s)
d = 0.18; % Separation between WFS array loudspeakers

% Noise source coefficient
amplitude = 1;
phase = 0;

% Frequencies
freqs = 0; % Random frequency. Required for technical purposes
numFreqs = 1; % Required for technical reasons

% Room characteristics and impulse response of chamber
WFS_AcPath_previously_calculated = false;
NS_AcPath_previously_calculated = true;
appendFreeSpaceAcPaths = false;
predefNumSampIR = false; % It will be calculated automatically
numReverbTime = 6;
beta = linspace(0, 0.5, numReverbTime); % Average reflection coefficient of the walls of the chamber

% WFS options
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

frequencyCorrection = true;

attenuationType = 'Miguel';

% Simulation options
timeDomainActive = true;
fakeTimeProcessing = false;
frequencyDomainActive = false;
automaticLengthModification = false;
predefSignals = true;
saveSignals = true;

SetupParametersScript
AcousticPathCalculationScript

simulationScript;

% Analysis
numSamp = size(rec_signals, 2);
freqsFFT = (0:numSamp - 1)*fs/numSamp;
sel = freqsFFT >= 0 & freqsFFT <= 1000;
freqsFFTsel = freqsFFT(sel);

axCanc = gobjects(numReverbTime, 1);
for rt = 1:numReverbTime
    recNS_freq = fft(recNS_signals(:, :, :, :, rt), [], 2); recNS_freq = recNS_freq(:, sel, :, :);
    rec_freq = fft(rec_signals(:, :, :, :, rt), [], 2); rec_freq = rec_freq(:, sel, :, :);
    canc = 20*log10(abs(rec_freq./recNS_freq));
    axCanc(rt) = histogram2D( canc, 2, freqsFFTsel, [], [] );
end

for rt = 1:numReverbTime
    axCanc(rt).XLabel.String = 'Frequency (Hz)';
    axCanc(rt).YLabel.String = 'Cancellation (dB)';
    axCanc(rt).YLim = [-20, 0];
    
    aux = num2str(beta(rt));
    if ismember('.', aux)
        aux = aux(3:end);
    end
    printfig(axCanc(rt).Parent, imagesPath, ['Experiment9_cancFreqCorrReflecCoef', aux], 'eps')
end

