%% Experiment 9
% Comparison of the WFS attenuations used by Miguel and by me.

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
f = 200;
NSsignal = sin(2*pi*f*t);

[POT_ad,POT_ar]=generamapa(-NSsignal', h_array, filtros_array, activo_array, NSsignal', h_ad_sources(:,1,:,:));

    %% B) Simulate the same scenario but with my code
roomDim = L;
predefRoomDim = true;
[X, Y] = ndgrid(mallado_x, mallado_y);
recPositions = [X(:), Y(:), 1.65*ones(numel(X), 1)];

WFSposition = [WFSpos', zPos*ones(96, 1)];
obj.WFSposition = WFSposition;
obj.WFSToolObj.scenarioObj.loudspeakersOrientation
broadsideDir = [cosd(tecta' - 90), sind(tecta' - 90), zeros(96, 1)];
obj.WFSToolObj.WFSarrayOrientation = simulator.vec2rotVec(broadsideDir);

WFSfilterLength = size(filtros_array, 1);
predefWFSfilterLength = true;

WFSarrayOffset = [0, 0, 0]; % [x, y, z] coordinates. Useful for generating acoustic path IR.

SetupParametersScript;
AcousticPathCalculationScript;

simulationScript;

%% Analysis of results

% Visualize

% Perform FFT
NS = fft(x');

    % My code
recWFS_signals = rec_signals - recNS_signals;
% recWFS = fft(recWFS_signals'); %/fs;
% rec = fft(rec_signals'); %/fs;
% recNS = fft(recNS_signals'); %/fs;
% H_WFS = recWFS./repmat(NS, [1, numMicro]);
% H_NS = recNS./repmat(NS, [1, numMicro]);
% corrFact = -recNS./recWFS;

    % Miguel's code
recWFS_signals_Miguel = [POT_ar', zeros(numMicro, length(x) - size(POT_ar, 1))]; % Mucho cuidado con esto!!! Al principio has interpretado que POT_ar era el campo recibido, pero es el campo proveniente del WFS!! Esto es importante. Te ha hecho perder varias horas
recNS_signals_Miguel = [POT_ad', zeros(numMicro, length(x) - size(POT_ad, 1))];
% recWFS_Miguel = fft(recWFS_signals_Miguel'); %/fs;
% rec_Miguel = fft(rec_signals_Miguel'); %/fs;
% recNS_Miguel = fft(recNS_signals_Miguel'); %/fs;
% H_WFS_Miguel = recWFS_Miguel./repmat(NS, [1, numMicro]);
% H_NS_Miguel = recNS_Miguel./repmat(NS, [1, numMicro]);
% corrFact_Miguel = -recNS_Miguel./recWFS_Miguel;

% Visualize
numSamp = size(recNS_signals, 2);
t = (0:numSamp-1)/fs;

axNS = axes(figure);
plot(axNS, t, recNS_signals(1,:), t, recNS_signals_Miguel(1,:))
axNS.XLabel.String = 'Time (s)';
axNS.YLabel.String = 'Signal (arbitrary units)';
axNS.Title.String = 'Received signal from the noise source';

axWFS = axes(figure);
plot(axWFS, t, recWFS_signals(1,:), t, d*recWFS_signals_Miguel(1,:))
axWFS.XLabel.String = 'Time (s)';
axWFS.YLabel.String = 'Signal (arbitrary units)';
axWFS.Title.String = 'Received signal from the WFS array';


% Why is the received signal from WFS different? The one received from the NS
% is the same!
wfsInd = 1;
recInd = 1;
[xInd, yInd] = ind2sub([length(mallado_x), length(mallado_y)], recInd);

tAcPath = (0:numSampIR - 1)/fs;
ax = axes(figure);
plot(ax, tAcPath, h_array(:, wfsInd, xInd, yInd), tAcPath, permute(WFSacPathIR(recInd, wfsInd, :), [3 1 2]))
plot(ax, tAcPath, h_array(:, wfsInd, xInd, yInd), tAcPath, permute(WFS_IR(recInd, :, wfsInd), [2 1 3]))

h_array_reshape = mergeAndPermute( h_array, {[3 4], 2 1});
isequal(WFSacPathIR, h_array_reshape) 
% Los caminos acústicos son iguales!! Muy bien :D

% Compare WFS filters
tWFSfilt = (0:WFSfilterLength - 1)/fs;

ax = axes(figure);
wfsInd = 2;
plot(ax, tWFSfilt, filtros_array(:, wfsInd), tWFSfilt, permute(obj.WFSToolObj.filtersWFS_IR(wfsInd, 1, :), [3 1 2]))
% There are differences of one sample in some deltas

% Let's test if that's the only remaining difference
obj.WFSToolObj.filtersWFS_IR = repmat(permute(filtros_array, [2 3 1]), [1 2 1]);
isequal(repmat(permute(filtros_array, [2 3 1]), [1 2 1]), obj.WFSToolObj.filtersWFS_IR)
obj.WFSToolObj.WFScalculation();
obj.WFSToolObj.simulate();
rec_signal2 = obj.WFSToolObj.simulField;
recWFS_signals2 = rec_signal2 - recNS_signals;

ax = axes(figure);
plot(ax, t, d*recWFS_signals2(1,:), t, recWFS_signals(1,:))
% La diferencia es muy poca, como era de esperar

axWFS = axes(figure);
plot(axWFS, t, recWFS_signals2(1,:), t, recWFS_signals_Miguel(1,:))
axWFS.XLabel.String = 'Time (s)';
axWFS.YLabel.String = 'Signal (arbitrary units)';
axWFS.Title.String = 'Received signal from the WFS array';

% Siguen sin ser iguales!!!

% [obj.WFSToolObj.scenarioObj.attenuations, an'] % The attenuations are the same

% Extraigo parte del código de generamapa.m
xsignal = -NSsignal';
POT_ar=zeros(length(NSsignal), numMicro);
WFSsignals = zeros(length(NSsignal), 96);
nx = size(h_array, 3);
ny = size(h_array, 4);
for altavoz=1:length(activo_array)
    cont=1;
    xfsignal=filter(filtros_array(:,activo_array(altavoz)),1,xsignal);
    WFSsignals(:, activo_array(altavoz)) = xfsignal;
    for yy=1:ny
        for xx=1:nx
        POT_ar(:,cont)=POT_ar(:,cont)+filter(h_array(:,activo_array(altavoz),xx,yy),1,xfsignal);
        cont=cont+1;
        end
    end
end
WFSsignals_ = [WFSsignals; zeros(length(x) - size(POT_ar, 1), 96)]';
size(WFSsignals_)
size(obj.WFSToolObj.WFSarrayCoefficient)

ax = axes(figure);
wfsInd = 96;
plot(ax, t, WFSsignals_(wfsInd, :), t, obj.WFSToolObj.WFSarrayCoefficient(wfsInd, :))

differ = WFSsignals' - obj.WFSToolObj.WFSarrayCoefficient(:, 1:44100);
max(abs(differ(:)))
% Las señales reproducidas por el WFS array son iguales!! ¿Por qué es el
% resultado diferente entonces?

obj.WFSToolObj.real = [false; false];
obj.WFSToolObj.WFScalculation();
obj.WFSToolObj.simulate();
recWFS_signal_special = obj.WFSToolObj.simulField;

axWFS = axes(figure);
plot(axWFS, t, recWFS_signal_special(1,:), t, recWFS_signals_Miguel(1,:))
axWFS.XLabel.String = 'Time (s)';
axWFS.YLabel.String = 'Signal (arbitrary units)';
axWFS.Title.String = 'Received signal from the WFS array';

max(abs(recWFS_signal_special(:) - recWFS_signals(:)))
% Son casi iguales, como debería de ser.
% Sigo sin entender por qué el resultado es diferente

% Voy a calcular la contribución de un solo altavoz a un solo micrófono
% Extraigo parte del código de generamapa.m
xsignal = -NSsignal';
POT_ar_red = zeros(length(NSsignal), 96);
nx = size(h_array, 3);
ny = size(h_array, 4);
for altavoz=1:length(activo_array)
    xfsignal = filter(filtros_array(:,activo_array(altavoz)), 1, xsignal);
    % max(abs(xfsignal' - obj.WFSToolObj.WFSarrayCoefficient(1, 1:44100)))

    POT_ar_red(:, activo_array(altavoz)) = filter(h_array(:,activo_array(altavoz),1,1),1,xfsignal);
end

recSignals_red = fftfilt_modified(permute(WFSacPathIR(1, :, :), [3 2 1]), obj.WFSToolObj.WFSarrayCoefficient')';


ax = axes(figure);
plot(ax, t, recSignals_red, t, [sum(POT_ar_red, 2); zeros(length(t) - length(NSsignal), 1)])
% Sí que da lo mismo así!!!!!!! No lo entiendo. ¿Cuál es la diferencia
% respecto al bucle original?
xsignal = -NSsignal';
POT_ar2 = zeros(length(NSsignal), numMicro);
nx = size(h_array, 3);
ny = size(h_array, 4);
for altavoz=1:length(activo_array)
    cont=1;
    xfsignal=filter(filtros_array(:,activo_array(altavoz)),1,xsignal);
    for yy=1:ny
        for xx=1:nx
            POT_ar2(:,cont)=POT_ar2(:,cont)+filter(h_array(:,activo_array(altavoz),xx,yy),1,xfsignal);
            cont=cont+1;
        end
    end
end
POT_ar2_ = [POT_ar2; zeros(length(t) - length(xsignal), numMicro)];

isequal(POT_ar2, POT_ar)

recSignals_lowLevel = fftfilt_modified(permute(WFSacPathIR, [3 2 1]), obj.WFSToolObj.WFSarrayCoefficient')';

% s = whos;
% bytes = [s.bytes];
% names = {s.name};
% [bytesSort, ind] = sort(bytes, 'descend');
% bytesSort(1:10)

ax = axes(figure);
plot(ax, t, recSignals_lowLevel, t, [sum(POT_ar_red, 2); zeros(length(t) - length(NSsignal), 1)])
plot(ax, t(1:44100), recWFS_signals2(1,1:44100), t(1:44100), POT_ar(:, 1))
plot(axWFS, t, recWFS_signals2(1,:), t, recWFS_signals_Miguel(1,:))


