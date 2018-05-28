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

%% Check that Miguel's scenario parameters produce the same result in both simulations

    %% A) Simulate with Miguel code
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

% Representamos el mapa de presiones
 dv=50; % Duración de la ventana para promediar al potencia
 despv=25;  % Desplazamiento del enventanado para realizar el calculo de la
% potencia

dibujapot(POT_ad+POT_ar ,L, WFSpos, NSpositions, mallado_x, mallado_y,dv,despv);



    %% B) Simulate the same scenario but with my code
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
