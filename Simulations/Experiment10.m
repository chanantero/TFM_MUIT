%% Experiment 10

% This file will try to solve the questions asked by the professor Miguel
% in the control session of 1st of June of 2018.

% - Are the acoustic path IR correct? What happens if you try different
% sampling frequencies? Will the temporal resolution be a problem?
% 

%% Preamble
pathSetUp;

imagesPath = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\TFM\Img\';

dataPathName = [globalPath, 'Data\'];
ID = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

%% Global parameters
% Precalculation of useful parameters

% Generate WFS array positions
s = WFSToolSimple.generateScenario(96);
WFSposDefault = s.loudspeakersPosition;
rotMat = [0 1 0; -1 0 0; 0 0 1];
WFSpositionRot = (rotMat * WFSposDefault')'; % Rotate 90º to the right
zPos = 1.65;
WFSoffset = [2.21 4.02 zPos];
WFSposition = WFSpositionRot + repmat(WFSoffset, 96, 1); % Add offset
broadsideDir = s.loudspeakersOrientation;
broadsideDir = (rotMat * broadsideDir')';

centreWFSx = (min(WFSposition(:, 1)) + max(WFSposition(:, 1)))/2;
centreWFSy = (min(WFSposition(:, 2)) + max(WFSposition(:, 2)))/2;

% Constants
roomDim = [9.13, 4.48, 2.64];
c = 340;


%% Section 1
% The questions we are trying to answer are:
% Are the acoustic path IR correct? What happens if you try different
% sampling frequencies? Will the temporal resolution be a problem?

% The acoustic paths are generated with the function rir_generator.

% A)
% Phase shift from theoretical to real case when the distance between
% source and receiver is fixed (one source and one receiver), for different
% sampling frequencies and different signal frequencies (x axis).

% Chooose a noise source position and a receiver position and a vector of
% sampling frequencies
NSposition = [1 1 zPos];
recPosition = [centreWFSx centreWFSy zPos];
fsVec = [4000, 8000, 16000, 32000, 44100];
numFS = length(fsVec);

dist = calcDistances(NSposition, recPosition);

% Generate impulse responses
acPathIR = cell(numFS, 1);
for k = 1:numFS
    fs = fsVec(k);
    numSampIR = floor(dist/c*fs * 2) + 1;
    acPathIR{k} = rir_generator(c, fs, recPosition, NSposition, roomDim, 0, numSampIR);
end

delayTheo = dist/c;
delays = zeros(numFS, 1);
phNumPoints = 44100;
phaseResponse = zeros(numFS, phNumPoints);
% fftNumPoints = 44100;
% acPathFR = cell(numFS, 1);
% phaseDelay = zeros(numFS, phdNumPoints);
% ax = axes(figure, 'NextPlot', 'Add');
for k = 1:numFS
    IR = acPathIR{k};
    numSamp = length(IR);
    fs = fsVec(k);
    tVec = (0:numSamp - 1)/fs;
    
    [~, ind] = max(IR);
    delays(k) = tVec(ind);
    
    
    [phi, ~] = phasez(IR, 1, phNumPoints);
    phaseResponse(k, :) = phi';
    
%     plot(ax, tVec, acPathIR{k})
    
%     acPathFR{k} = fft(IR, fftNumPoints);
%     [phi, w] = phasedelay(IR, 1, phdNumPoints);
%     phaseDelay(k, :) = phi'/fs;
end
% delays - delayTheo

% % Phase delay response: 2 ways of doing it with different result. Why??
% % What does the phasedelay function do exactly?
% ax = axes(figure, 'NextPlot', 'Add');
% for k = 1:numFS
%     fs = fsVec(k);
%     f = w/(2*pi)*fs;
%     plot(ax, f, phaseDelay(k, :))
% end
% hold on
% plot(ax, f, delayTheo*ones(1, length(f)))
% ax.XLim = [0 1000];
% ax.YLim = [0 dist/c*2];
% ax.XLabel.String = 'Frequency (Hz)';
% 
ax = axes(figure, 'NextPlot', 'Add');
for k = 1:numFS
    fs = fsVec(k);
    f = w/(2*pi)*fs;
    phaseDelayCustom = -phaseResponse(k,:)./(2*pi*f');
    plot(ax, f, phaseDelayCustom)
end
plot(ax, f, delayTheo*ones(1, length(f)))
ax.XLim = [0 1000];
ax.YLim = [0 dist/c*2];
ax.XLabel.String = 'Frequency (Hz)';

ax = axes(figure, 'NextPlot', 'Add');
for k = 1:numFS
    fs = fsVec(k);
    f = w/(2*pi)*fs;
    idealPhaseResponse = -f*dist/c*2*pi;
    plot(ax, f, rad2deg(phaseResponse(k, :) - idealPhaseResponse'))    
end
ax.XLim = [0 1000];
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = '\angle{H_{real}} - \angle{H_{ideal}} (º)';

% Filter a tone to see the effects of the filters
durSign = 1;
f0 = 800;
filtX = cell(numFS, 1);
ax = axes(figure, 'NextPlot', 'Add');
for k = 1:numFS
    fs = fsVec(k);
    N = floor(durSign*fs);
    t = (0:N - 1)/fs;
    x = cos(2*pi*f0*t);
    
     x_f = filter(acPathIR{k}, 1, x);
     
     filtX{k} = x_f;
     
     plot(ax, t, x_f)
end
filtXideal = cos(2*pi*f0*(t - delayTheo));
plot(ax, t, filtXideal*0.02, '--');

% Conclussion: rir_generator takes in account the temporal resolution and
% compensates for it. The phase response is basically the same for all
% sampling frequencies, although it is shifted from the ideal case
% depending on the signal frequency.

% B) The same that A) but with different positions of the receiver.
% The aim is to see whether this signal-frequency-dependent phase-shift is 
% constant with the distance

% Chooose a noise source position and different receiver positions and 
% sampling frequencies
NSposition = [1 1 zPos];
numRec = 10;
recPosX = centreWFSx + linspace(-3, 6, numRec)';
recPosition = [recPosX repmat([centreWFSy zPos], numRec, 1)];
fsVec = [4000, 8000, 16000, 32000, 44100]';
numFS = length(fsVec);

dist = calcDistances(NSposition, recPosition);

% Generate impulse responses
acPathIR = cell(numFS, 1);
for k = 1:numFS
    fs = fsVec(k);
    numSampIR = floor(max(dist)/c*fs * 2) + 1;
    acPathIR{k} = rir_generator(c, fs, recPosition, NSposition, roomDim, 0, numSampIR);
end

delayTheo = dist/c;
phNumPoints = 44100;
phaseResponse = cell(numFS, 1);
for k = 1:numFS
    IRs = acPathIR{k};
    numSamp = size(IRs, 2);
    fs = fsVec(k);
    tVec = (0:numSamp - 1)/fs;
    phaseResponseCurr = zeros(numRec, phNumPoints);
    for r = 1:numRec
        phi = phasez(IRs(r,:), 1, phNumPoints);
        phaseResponseCurr(r, :) = phi';
    end
    phaseResponse{k} = phaseResponseCurr;
end
fVecs = repmat((0:phNumPoints - 1)/(phNumPoints*2), numFS, 1).*...
    repmat(fsVec, 1, phNumPoints);

ax = axes(figure, 'NextPlot', 'Add');
lin = gobjects(numFS, numRec);
for k = 1:numFS
    f = fVecs(k, :);
    phaseResponseCur = phaseResponse{k};
    for r = 1:numRec
        idealPhaseResponse = -f*dist(r)/c*2*pi;
        lin(k, r) = plot(ax, f, rad2deg(phaseResponseCur(r, :) - idealPhaseResponse));
    end
end
ax.XLim = [0 1000];
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = '\angle{H_{real}} - \angle{H_{ideal}} (º)';

activeFS = 1:5;
aux = repmat({'off'}, size(lin));
[lin(:).Visible] = aux{:};
aux = repmat({'on'}, numRec*numel(activeFS), 1);
[lin(activeFS, :).Visible] = aux{:};

activeRec = [7];
aux = repmat({'off'}, size(lin));
[lin(:).Visible] = aux{:};
aux = repmat({'on'}, numFS*numel(activeRec), 1);
[lin(:, activeRec).Visible] = aux{:};

% Conclussion: the phase shift is constant with sampling frequency and distance and
% so, it doesn't affect the final calculation of cancellation.

%% Section 2: Compare acoutic path IR in near field with far field
% The purpose is to see if rir_generator performs some kind of near far
% field approximation, because maybe the euse of the far field condition is
% not realistic

zPos = 1.65;
roomDim = [9.13, 4.48, 2.64];
c = 340;

% Chooose a noise source position and different receiver positions and 
% sampling frequencies
NSposition = [1 1 zPos];
numRec = 100;
recPosX = NSposition(1) + linspace(0.01, 20, numRec)';
recPosition = [recPosX repmat([NSposition(2) zPos], numRec, 1)];
fs = 44100;

dist = calcDistances(NSposition, recPosition);

% Generate impulse responses
numSampIR = floor(max(dist)/c*fs * 2) + 1;
acPathIR = rir_generator(c, fs, recPosition, NSposition, roomDim, 0, numSampIR);

% Get amplitude of filter depending on frequency
acPathFR = fft(acPathIR');
fVec = (0:numSampIR - 1)*fs/numSampIR;

% sel = fVec >= 0 & fVec <= 1000;
% acPathFRsel = acPathFR(sel, :);
% visualObj = animation({fVec(sel), 1:numRec},...
%     {abs(acPathFRsel)}, {'Frequency', 'Receiver point'}, ...
%     {'Frequency response'}, [], []);

selFreq = [100 200 400 600 800];
acPathSel = interp1(fVec, acPathFR, selFreq);

% 1st way of curve fitting
% data.x = dist';
% a = zeros(numel(selFreq), 1);
% for f = 1:numel(selFreq)
%     data.y = abs(acPathSel(f, :))';
%     [ params, gofs ] = fitInterface( data, {'a/x'} );
%     a(f) = params{1}.a;
% end

% 2nd way of curve fitting
normTheoLog = log(1./dist);
data.x = dist';
aLog = zeros(numel(selFreq), 1);
for f = 1:numel(selFreq)
    data.y = log(abs(acPathSel(f, :))');
    [ params, gofs ] = fitInterface( data, {'a - log(x)'} );
    aLog(f) = params{1}.a;
end
a = exp(aLog);

[aExt, distExt] = pointWiseExtend(a', dist');

ax = axes(figure, 'NextPlot', 'Add');
lFR = plot(ax, dist, abs(acPathSel)');
lTheoFit = plot(ax, dist, aExt./distExt);

indFreq = 5;
aux = repmat({'off'}, size(lFR));
[lFR(:).Visible] = aux{:};
[lTheoFit(:).Visible] = aux{:};
aux = repmat({'on'}, numel(indFreq), 1);
[lFR(indFreq).Visible] = aux{:};
[lTheoFit(indFreq).Visible] = aux{:};

absTheo = aExt./distExt;
rel = abs(acPathSel)'./absTheo;
relLog = log(rel);
exp(relLog(:, 1))

ax = axes(figure, 'NextPlot', 'Add');
plot(ax, dist, relLog);

% Plot IR of the closest receiver point
ax = axes(figure);
plot(ax, acPathIR(1,:))

% Conclussion: although it is not exact, the acoustic paths IR generated by
% rir_generator follow approximatelly the a/dist curve form.

%% Section 3: How does the sampling frequency influence performance?
% As it has been proved in section 1, the acoustic paths IR generated by
% rir_generator seem to have a very good response even for low sampling
% frequencies.
% Then, the source of error comes necessarily from the generation of WFS
% signals.
% The way we are going to evaluate this error is the next.
% For each WFS element, a signal is created by the filtering of the noise
% source signal. As we can calculate exactly the amount of delay produced
% by each WFS filter, it can be compared with the theoretical one.
% With that information we can calculate de phase shift for each frequency
zPos = 1.65;
roomDim = [9.13, 4.48, 2.64];
c = 340;

% Chooose a noise source position and different receiver positions and 
% sampling frequencies
NSposition = [1 1 zPos];
numRec = 100;
recPosX = NSposition(1) + linspace(0.01, 20, numRec)';
recPosition = [recPosX repmat([NSposition(2) zPos], numRec, 1)];
fs = 44100;

dist = calcDistances(NSposition, recPosition);

delaysTheo = dist/c;

indDelta = floor(delaysTheo*fs) + 1;
delaysReal = (indDelta - 1)/fs;

max(abs(delaysReal - delaysTheo)) <= 1/fs


%% Section 4: Far field behavior
% I expect that the correction factor when the noise sources are far from
% the WFS array don't change too much with respect when they are close.
% However, this test will clarify it.

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
predefWFSfilterLength = true;
% WFSfilterLength is calculated some lines below

% Frequencies
freqs = [20:20:1000];
% freqs = 0;

% Room characteristics and impulse response of chamber
WFS_AcPath_previously_calculated = true;
NS_AcPath_previously_calculated = true;
appendFreeSpaceAcPaths = true;
predefNumSampIR = false;
beta = [];
WFSarrayOffset = [0 0 0];    

% WFS options
frequencyCorrection = false;
attenuationType = 'Miguel';

% Simulation options
timeDomainActive = false;
fakeTimeProcessing = false;
frequencyDomainActive = true;
automaticLengthModification = false;
predefSignals = true;
saveSignals = false;

durSign = 1; % Duration of signal
t = (0:ceil(durSign*fs)-1)/fs;
NSsignal = chirp(t, 20, durSign, 940);

% Positions of the noise source
% Quarter of a circle
numPointsPerArc = 10;
radius = [100];
numArcs = numel(radius);
centreX = (max(WFSposition(:, 1)) + min(WFSposition(:, 1)))/2;
centreY = (max(WFSposition(:, 2)) + min(WFSposition(:, 2)))/2;
y1 = centreY;
y2 = roomDim(2) - centreY;
alphaMax = pi + asin(y1/max(radius)) - deg2rad(1);
alphaMin = pi - asin(y2/max(radius)) + deg2rad(1);
alpha = linspace(alphaMin, alphaMax, numPointsPerArc)';
x = centreX + repmat(radius, numPointsPerArc, 1).*repmat(cos(alpha), 1, numArcs);
y = centreY + repmat(radius, numPointsPerArc, 1).*repmat(sin(alpha), 1, numArcs);
NSpositions = [x(:), y(:), zPos*ones(numel(x), 1)];

% Microphone positions
% Rectangular grid
marginRatio = 0.3;
numPointsX = 5;
numPoinstY = 3;
extRectXmin = min(WFSposition(:, 1));
extRectXmax = max(WFSposition(:, 1));
extRectYmin = min(WFSposition(:, 2));
extRectYmax = max(WFSposition(:, 2));
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
z = zPos;
[X, Y, Z] = ndgrid(x, y, z);
recPositions = [X(:), Y(:), Z(:)];

predefRoomDim = true;

if ~exist('obj', 'var') || ~isvalid(obj)
    obj = SimulationController;
end

obj.WFSposition = WFSposition;
obj.WFSToolObj.WFSarrayOrientation = simulator.vec2rotVec(broadsideDir);
obj.WFSToolObj.scenarioObj.ax.XLim = [0, roomDim(1)];
obj.WFSToolObj.scenarioObj.ax.YLim = [0, roomDim(2)];

dist = calcDistances(WFSposition, NSpositions);
WFSfilterLength = (floor(max(dist(:))/c*fs) + 1)*2;

SetupParametersScript
AcousticPathCalculationScript
simulationScript;

% Calculate individual and global correction factors
sFreq = s(2, :, :, :);
[sFreq, corrFactInd, corrFactGlobal] = SimulationController.addCancellationParametersToStructure(sFreq);

% Hisgoram
absCorrFactEdges = 0:0.1:4;
phaseCorrFactEdges = -180:180;

axAbsInd = histogram2D( abs(corrFactInd), 2, freqs, [], absCorrFactEdges );
axPhaseInd = histogram2D( rad2deg(angle(corrFactInd)), 2, freqs, [], [] );
% axPhaseInd.YLim = [0 50];

axAbsGlob = histogram2D( abs(corrFactGlobal), 2, freqs, [], [] );
axPhaseGlob= histogram2D( rad2deg(angle(corrFactGlobal)), 2, freqs, [], []);
% axPhaseGlob.YLim = [0 50];

axAbsInd.Parent.Name = 'CorrFact Ind Abs';
axPhaseInd.Parent.Name = 'CorrFact Ind Phase';
axAbsGlob.Parent.Name = 'CorrFact Glob Abs';
axPhaseGlob.Parent.Name = 'CorrFact Glob Phase';

corrFactTheo = sqrt(1i*freqs/c);

axAbsInd.NextPlot = 'Add';
plot(axAbsInd, freqs, abs(corrFactTheo), 'r', 'LineWidth', 3);
axAbsGlob.NextPlot = 'Add';
plot(axAbsGlob, freqs, abs(corrFactTheo), 'r', 'LineWidth', 3);

% Conclussion: the tendency is the same in far-field

%% Section 5: Why does the phase shift of the correction factor change linearly for low frequencies?
% This is pretty misterious.
% My first guess is that the WFS equations, since they are deduced from
% kirchhoff's equation when the distance between noise source and the
% cancellation volume is higher than some number of wavelengths, don't work
% for higher frequencies.

% Use configuration in Section 4.

% Normal configuration
d = 0.18; % Separation between two contiguous loudspeakers. Size of one loudspeaker
nb = 8; % Bottom and upper sides of the octogon (2 sides)
nd = 8; % Diagonal sides of the octogon (4 sides)
nl = 24; % Lateral side of the octogon (2 sides)
betabd = 45; % Deviation angle between bottom/upper and diagonal sides

[ x, y, alfa ] = octogon(d, nb, nd, nl, betabd);
z = zeros(numel(x), 1);
WFSposDefault = [x, y, z];
broadsideDir = [cosd(alfa), sind(alfa), zeros(numel(alfa), 1)];
                    
rotMat = [0 1 0; -1 0 0; 0 0 1];
WFSpositionRot = (rotMat * WFSposDefault')'; % Rotate 90º to the right
zPos = 1.65;
WFSoffset = [2.21 4.02 zPos];
WFSposition = WFSpositionRot + repmat(WFSoffset, 96, 1); % Add offset
broadsideDir = (rotMat * broadsideDir')';

obj.WFSposition = WFSposition;
obj.WFSToolObj.WFSarrayOrientation = simulator.vec2rotVec(broadsideDir);
obj.WFSToolObj.scenarioObj.ax.XLim = [0, roomDim(1)];
obj.WFSToolObj.scenarioObj.ax.YLim = [0, roomDim(2)];
obj.WFSToolObj.WFSarrayAdjacentSeparation = d;

% Microphone positions
% Rectangular grid
marginRatio = 0.2;
numPointsX = 5;
numPoinstY = 3;
extRectXmin = min(WFSposition(:, 1));
extRectXmax = max(WFSposition(:, 1));
extRectYmin = min(WFSposition(:, 2));
extRectYmax = max(WFSposition(:, 2));
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
z = zPos;
[X, Y, Z] = ndgrid(x, y, z);
recPositions = [X(:), Y(:), Z(:)];

% Positions of the noise source
% Quarter of a circle
numPointsPerArc = 2;
radius = [100];
numArcs = numel(radius);
centreX = (max(WFSposition(:, 1)) + min(WFSposition(:, 1)))/2;
centreY = (max(WFSposition(:, 2)) + min(WFSposition(:, 2)))/2;
y1 = centreY;
y2 = roomDim(2) - centreY;
alphaMax = pi/2;
alphaMin = 0;
alpha = linspace(alphaMin, alphaMax, numPointsPerArc)';
x = centreX + repmat(radius, numPointsPerArc, 1).*repmat(cos(alpha), 1, numArcs);
y = centreY + repmat(radius, numPointsPerArc, 1).*repmat(sin(alpha), 1, numArcs);
NSpositions = [x(:), y(:), zPos*ones(numel(x), 1)];

SetupParametersScript
AcousticPathCalculationScript
simulationScript;

% Calculate individual and global correction factors
sFreq = s(2, :, :, :);
[sFreq, corrFactInd, corrFactGlobal] = SimulationController.addCancellationParametersToStructure(sFreq);

% Hisgoram
absCorrFactEdges = 0:0.1:4;
phaseCorrFactEdges = -180:180;

axAbsInd = histogram2D( abs(corrFactInd), 2, freqs, [], absCorrFactEdges );
axPhaseInd = histogram2D( rad2deg(angle(corrFactInd)), 2, freqs, [], [] );
% axPhaseInd.YLim = [0 50];

axAbsGlob = histogram2D( abs(corrFactGlobal), 2, freqs, [], [] );
axPhaseGlob= histogram2D( rad2deg(angle(corrFactGlobal)), 2, freqs, [], []);
% axPhaseGlob.YLim = [0 50];

axAbsInd.Parent.Name = 'CorrFact Ind Abs';
axPhaseInd.Parent.Name = 'CorrFact Ind Phase';
axAbsGlob.Parent.Name = 'CorrFact Glob Abs';
axPhaseGlob.Parent.Name = 'CorrFact Glob Phase';

corrFactTheo = sqrt(1i*freqs/c);

axAbsInd.NextPlot = 'Add';
plot(axAbsInd, freqs, abs(corrFactTheo), 'r', 'LineWidth', 3);
axAbsGlob.NextPlot = 'Add';
plot(axAbsGlob, freqs, abs(corrFactTheo), 'r', 'LineWidth', 3);


% B) Let's use a WFS array with a separation between loudspeakers d = 0.9
d = 0.36; % Separation between two contiguous loudspeakers. Size of one loudspeaker
nb = 8; % Bottom and upper sides of the octogon (2 sides)
nd = 8; % Diagonal sides of the octogon (4 sides)
nl = 24; % Lateral side of the octogon (2 sides)
betabd = 45; % Deviation angle between bottom/upper and diagonal sides

[ x, y, alfa ] = octogon(d, nb, nd, nl, betabd);
z = zeros(numel(x), 1);
WFSposDefault = [x, y, z];
broadsideDir = [cosd(alfa), sind(alfa), zeros(numel(alfa), 1)];
                    
rotMat = [0 1 0; -1 0 0; 0 0 1];
WFSpositionRot = (rotMat * WFSposDefault')'; % Rotate 90º to the right
zPos = 1.65;
WFSoffset = [2.21 4.02 zPos];
WFSposition = WFSpositionRot + repmat(WFSoffset, 96, 1); % Add offset
broadsideDir = (rotMat * broadsideDir')';

obj.WFSposition = WFSposition;
obj.WFSToolObj.WFSarrayOrientation = simulator.vec2rotVec(broadsideDir);
obj.WFSToolObj.scenarioObj.ax.XLim = [0, roomDim(1)];
obj.WFSToolObj.scenarioObj.ax.YLim = [0, roomDim(2)];
obj.WFSToolObj.WFSarrayAdjacentSeparation = d;

% Microphone positions
% Rectangular grid
marginRatio = 0.2;
numPointsX = 5;
numPoinstY = 3;
extRectXmin = min(WFSposition(:, 1));
extRectXmax = max(WFSposition(:, 1));
extRectYmin = min(WFSposition(:, 2));
extRectYmax = max(WFSposition(:, 2));
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
z = zPos;
[X, Y, Z] = ndgrid(x, y, z);
recPositions = [X(:), Y(:), Z(:)];

% Positions of the noise source
% Quarter of a circle
numPointsPerArc = 2;
radius = [100];
numArcs = numel(radius);
centreX = (max(WFSposition(:, 1)) + min(WFSposition(:, 1)))/2;
centreY = (max(WFSposition(:, 2)) + min(WFSposition(:, 2)))/2;
y1 = centreY;
y2 = roomDim(2) - centreY;
alphaMax = pi/2;
alphaMin = 0;
alpha = linspace(alphaMin, alphaMax, numPointsPerArc)';
x = centreX + repmat(radius, numPointsPerArc, 1).*repmat(cos(alpha), 1, numArcs);
y = centreY + repmat(radius, numPointsPerArc, 1).*repmat(sin(alpha), 1, numArcs);
NSpositions = [x(:), y(:), zPos*ones(numel(x), 1)];

SetupParametersScript
AcousticPathCalculationScript
simulationScript;

% Calculate individual and global correction factors
sFreq = s(2, :, :, :);
[sFreq, corrFactInd, corrFactGlobal] = SimulationController.addCancellationParametersToStructure(sFreq);

% Hisgoram
absCorrFactEdges = 0:0.1:4;
phaseCorrFactEdges = -180:180;

axAbsInd = histogram2D( abs(corrFactInd), 2, freqs, [], absCorrFactEdges );
axPhaseInd = histogram2D( rad2deg(angle(corrFactInd)), 2, freqs, [], [] );
% axPhaseInd.YLim = [0 50];

axAbsGlob = histogram2D( abs(corrFactGlobal), 2, freqs, [], [] );
axPhaseGlob= histogram2D( rad2deg(angle(corrFactGlobal)), 2, freqs, [], []);
% axPhaseGlob.YLim = [0 50];

axAbsInd.Parent.Name = 'CorrFact Ind Abs';
axPhaseInd.Parent.Name = 'CorrFact Ind Phase';
axAbsGlob.Parent.Name = 'CorrFact Glob Abs';
axPhaseGlob.Parent.Name = 'CorrFact Glob Phase';

corrFactTheo = sqrt(1i*freqs/c);

axAbsInd.NextPlot = 'Add';
plot(axAbsInd, freqs, abs(corrFactTheo), 'r', 'LineWidth', 3);
axAbsGlob.NextPlot = 'Add';
plot(axAbsGlob, freqs, abs(corrFactTheo), 'r', 'LineWidth', 3);

% Interesting, the transition zone doubles when size of the enclosed area
% (with the loudspeaker separation) gets reduced to a half.

% Let's keep the original size, but double the number of loudspeakers
dI = 20; % Density increment
d = 0.18/dI; % Separation between two contiguous loudspeakers. Size of one loudspeaker
nb = 8*dI; % Bottom and upper sides of the octogon (2 sides)
nd = 8*dI; % Diagonal sides of the octogon (4 sides)
nl = 24*dI; % Lateral side of the octogon (2 sides)
betabd = 45; % Deviation angle between bottom/upper and diagonal sides

[ x, y, alfa ] = octogon(d, nb, nd, nl, betabd);
z = zeros(numel(x), 1);
WFSposDefault = [x, y, z];
broadsideDir = [cosd(alfa), sind(alfa), zeros(numel(alfa), 1)];
                    
rotMat = [0 1 0; -1 0 0; 0 0 1];
WFSpositionRot = (rotMat * WFSposDefault')'; % Rotate 90º to the right
zPos = 1.65;
WFSoffset = [2.21 4.02 zPos];
WFSposition = WFSpositionRot + repmat(WFSoffset, 96*dI, 1); % Add offset
broadsideDir = (rotMat * broadsideDir')';

obj.WFSposition = WFSposition;
obj.WFSToolObj.WFSarrayOrientation = simulator.vec2rotVec(broadsideDir);
obj.WFSToolObj.scenarioObj.ax.XLim = [0, roomDim(1)];
obj.WFSToolObj.scenarioObj.ax.YLim = [0, roomDim(2)];
obj.WFSToolObj.WFSarrayAdjacentSeparation = d;

% Microphone positions
% Rectangular grid
marginRatio = 0.2;
numPointsX = 5;
numPoinstY = 3;
extRectXmin = min(WFSposition(:, 1));
extRectXmax = max(WFSposition(:, 1));
extRectYmin = min(WFSposition(:, 2));
extRectYmax = max(WFSposition(:, 2));
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
z = zPos;
[X, Y, Z] = ndgrid(x, y, z);
recPositions = [X(:), Y(:), Z(:)];

% Positions of the noise source
% Quarter of a circle
numPointsPerArc = 2;
radius = [100];
numArcs = numel(radius);
centreX = (max(WFSposition(:, 1)) + min(WFSposition(:, 1)))/2;
centreY = (max(WFSposition(:, 2)) + min(WFSposition(:, 2)))/2;
y1 = centreY;
y2 = roomDim(2) - centreY;
alphaMax = pi/2;
alphaMin = 0;
alpha = linspace(alphaMin, alphaMax, numPointsPerArc)';
x = centreX + repmat(radius, numPointsPerArc, 1).*repmat(cos(alpha), 1, numArcs);
y = centreY + repmat(radius, numPointsPerArc, 1).*repmat(sin(alpha), 1, numArcs);
NSpositions = [x(:), y(:), zPos*ones(numel(x), 1)];

SetupParametersScript
AcousticPathCalculationScript
simulationScript;

% Calculate individual and global correction factors
sFreq = s(2, :, :, :);
[sFreq, corrFactInd, corrFactGlobal] = SimulationController.addCancellationParametersToStructure(sFreq);

% Hisgoram
absCorrFactEdges = 0:0.1:4;
phaseCorrFactEdges = -180:180;

axAbsInd = histogram2D( abs(corrFactInd), 2, freqs, [], [] );
axPhaseInd = histogram2D( rad2deg(angle(corrFactInd)), 2, freqs, [], [] );
% axPhaseInd.YLim = [0 50];

axAbsGlob = histogram2D( abs(corrFactGlobal), 2, freqs, [], [] );
axPhaseGlob= histogram2D( rad2deg(angle(corrFactGlobal)), 2, freqs, [], []);
% axPhaseGlob.YLim = [0 50];

axAbsInd.Parent.Name = 'CorrFact Ind Abs';
axPhaseInd.Parent.Name = 'CorrFact Ind Phase';
axAbsGlob.Parent.Name = 'CorrFact Glob Abs';
axPhaseGlob.Parent.Name = 'CorrFact Glob Phase';

corrFactTheo = sqrt(1i*freqs/c);

axAbsInd.NextPlot = 'Add';
plot(axAbsInd, freqs, abs(corrFactTheo), 'r', 'LineWidth', 3);
axAbsGlob.NextPlot = 'Add';
plot(axAbsGlob, freqs, abs(corrFactTheo), 'r', 'LineWidth', 3);

% Conclussion: the transition zone doesn't depend on the separation between
% loudspeakers

% It seems that the critical factor is the physical size in meters of the
% WFS array. Maybe it has to do with the far field condition. This is, the
% critical factor might be the distance in wavelengths from the receiving 
% point to the
% line of secondary sources.

% In order to go deep in tis concept, I'm going to create a truncated
% linear array, and I'm going to play with the longitude of it and the
% point where the receiver point is situated
numLongWFS = 100;
longWFS = 0.18*23*linspace(0.1, 10, numLongWFS);
zPos = 1.65;
numWFS = 200;
WFSyPos = 3;

numRec = 100;
distWFS2Rec = linspace(0.1, 10, numRec);

yNS = (10:10:100)';
numNS = length(yNS);

progressBarActive = false;

ss = cell(numel(longWFS), 1);
for long = 1:numel(longWFS)
    fprintf('long = %d/%d\n', long, numel(longWFS));
    
d = longWFS(long)/numWFS;
WFSposition = [(0:numWFS-1)'*d, WFSyPos*ones(numWFS, 1), ones(numWFS, 1)*zPos];
WFScentreX = (max(WFSposition(:, 1)) + min(WFSposition(:, 1)))/2;
broadsideDir = repmat([0 -1 0], numWFS, 1);

recPositions = [WFScentreX*ones(numRec, 1), WFSyPos - distWFS2Rec', zPos*ones(numRec, 1)];

NSpositions = [WFScentreX*ones(numNS, 1), yNS, zPos*ones(numNS, 1)];

obj.WFSposition = WFSposition;
obj.WFSToolObj.WFSarrayOrientation = simulator.vec2rotVec(broadsideDir);
obj.WFSToolObj.scenarioObj.ax.XLim = [min(WFSposition(:, 1)), max(WFSposition(:, 1))];
obj.WFSToolObj.scenarioObj.ax.YLim = [min(recPositions(:, 2)), NSpositions(2)];
obj.WFSToolObj.WFSarrayAdjacentSeparation = d;

SetupParametersScript
AcousticPathCalculationScript
simulationScript;

sFreq = s(2, :, :, :);
ss{long} = sFreq;
end

sFreqTotal = cat(1, ss{:});

% Calculate individual and global correction factors
[sFreq, corrFactInd, corrFactGlobal] = SimulationController.addCancellationParametersToStructure(sFreqTotal);

% Hisgoram
absCorrFactEdges = 0:0.1:4;
phaseCorrFactEdges = -180:180;

axAbsInd = histogram2D( abs(corrFactInd), 2, freqs, [], [] );
axPhaseInd = histogram2D( rad2deg(angle(corrFactInd)), 2, freqs, [], [] );
% axPhaseInd.YLim = [0 50];

axAbsGlob = histogram2D( abs(corrFactGlobal), 2, freqs, [], [] );
axPhaseGlob= histogram2D( rad2deg(angle(corrFactGlobal)), 2, freqs, [], []);
% axPhaseGlob.YLim = [0 50];

axAbsInd.Parent.Name = 'CorrFact Ind Abs';
axPhaseInd.Parent.Name = 'CorrFact Ind Phase';
axAbsGlob.Parent.Name = 'CorrFact Glob Abs';
axPhaseGlob.Parent.Name = 'CorrFact Glob Phase';

corrFactTheo = sqrt(1i*freqs/c);

axAbsInd.NextPlot = 'Add';
plot(axAbsInd, freqs, abs(corrFactTheo), 'r', 'LineWidth', 3);
axAbsGlob.NextPlot = 'Add';
plot(axAbsGlob, freqs, abs(corrFactTheo), 'r', 'LineWidth', 3);

ax = axes(figure);
plot(ax, squeeze(rad2deg(angle(corrFactInd))))

visualObj = animation({longWFS, freqs, 1:numRec},...
    {squeeze(rad2deg(angle(corrFactInd)))}, {'Longitude WFSarray', 'Frequency', 'Receiver point'}, ...
    {'Frequency response'}, [], []);


%% Theoretical comprobations

C = linspace(1, 20, 50); % Escalation factor
numC = numel(C);
longWFS0 = 0.18*100;
longWFS = longWFS0*C;
zPos = 1.65;
numWFS = 200;
WFSyPos = 3;
broadsideDir = repmat([0 -1 0], numWFS, 1);

distWFS2Rec0 = 1;
distWFS2Rec = distWFS2Rec0*C;

distNS0 = 5;
yNS = WFSyPos + distNS0*C';
numNS = length(yNS);

freq0 = 200;

progressBarActive = false;

freqs = freq0;
ss = cell(numC, 1);
for esc = 1:numC
    
d = longWFS(esc)/numWFS;
WFSposition = [(0:numWFS-1)'*d, WFSyPos*ones(numWFS, 1), ones(numWFS, 1)*zPos];
WFScentreX = (max(WFSposition(:, 1)) + min(WFSposition(:, 1)))/2;

recPositions = [WFScentreX, WFSyPos - distWFS2Rec(esc), zPos];

NSpositions = [WFScentreX, yNS(esc), zPos];

obj.WFSposition = WFSposition;
obj.WFSToolObj.WFSarrayOrientation = simulator.vec2rotVec(broadsideDir);
obj.WFSToolObj.scenarioObj.ax.XLim = [min(WFSposition(:, 1)), max(WFSposition(:, 1))];
obj.WFSToolObj.scenarioObj.ax.YLim = [min(recPositions(:, 2)), 1+NSpositions(2)];
obj.WFSToolObj.WFSarrayAdjacentSeparation = d;

SetupParametersScript
AcousticPathCalculationScript
simulationScript;

sFreq = s(2, :, :, :);
ss{esc} = sFreq;
end


freqEsc = freq0*C;

d = longWFS0/numWFS;
WFSposition = [(0:numWFS-1)'*d, WFSyPos*ones(numWFS, 1), ones(numWFS, 1)*zPos];
WFScentreX = (max(WFSposition(:, 1)) + min(WFSposition(:, 1)))/2;

recPositions = [WFScentreX, WFSyPos - distWFS2Rec0, zPos];

NSpositions = [WFScentreX, WFSyPos + distNS0, zPos];

obj.WFSposition = WFSposition;
obj.WFSToolObj.WFSarrayOrientation = simulator.vec2rotVec(broadsideDir);
obj.WFSToolObj.scenarioObj.ax.XLim = [min(WFSposition(:, 1)), max(WFSposition(:, 1))];
obj.WFSToolObj.scenarioObj.ax.YLim = [min(recPositions(:, 2)), 1+NSpositions(2)];
obj.WFSToolObj.WFSarrayAdjacentSeparation = d;

sOr = cell(numC, 1);
for esc = 1:numC
freqs = freqEsc(esc);

SetupParametersScript
AcousticPathCalculationScript
simulationScript;

sFreq = s(2, :, :, :);
sOr{esc} = sFreq;
end

sFreqTotal = cat(1, ss{:});
sFreqOrTotal = cat(1, sOr{:});

recWFScoef = zeros(size(sFreqTotal));
recWFScoef(:) = [sFreqTotal.recWFScoef];

recWFScoefOr = zeros(size(sFreqOrTotal));
recWFScoefOr(:) = [sFreqOrTotal.recWFScoef];

ax = axes(figure);
plot(ax, C, abs([recWFScoef, recWFScoefOr./C']));

% Conclussion: scaling by C is the same as having the original scenario,
% calculate the field when k is k_original*C, and dividing it by C

%% Continuation of Secion 5: Structured Experiments

% A) Ideal case: L=inf, Q_perfect and Q
% Monopole source at (0, d_ps, 0)
% Infinite linear array at the x axis

limx = [-1000 1000]; % Limits of integration on the infinite source

numDist_rec = 100;
d_rec = linspace(0, 10, numDist_rec + 1); d_rec = d_rec(2:end);
x_rec = 0;
d_ps = linspaceExtended([0 1 10 100], [10, 9, 3]); d_ps = d_ps(2:end);
k = 1;

numDist_rec = length(d_rec);
numDist_ps = length(d_ps);
numX_rec = length(x_rec);

r_ps = @(d_ps, x) sqrt(d_ps.^2 + x.^2); % Distance from primary source to secondary source
r_rec = @(d_rec, x_rec, x_s) sqrt(d_rec.^2 + (x_rec - x_s).^2); % Distance from point of the secondary source and point of measure
Q_perfect = @(x_s, d_ps) d_ps./r_ps(d_ps, x_s).*exp(-1i*k*r_ps(d_ps, x_s))./sqrt(r_ps(d_ps, x_s)).*(sqrt(1i*k/(2*pi)) + 1./(r_ps(d_ps, x_s)*sqrt(1i*2*pi*k))); % Feeding of each infinitesimal secondary source
Q = @(x_s, d_ps) d_ps./r_ps(d_ps, x_s).*exp(-1i*k*r_ps(d_ps, x_s))./sqrt(r_ps(d_ps, x_s)).*sqrt(1i*k/(2*pi)); % Feeding of each infinitesimal secondary source
G = @(dist, k) exp(-1i*k*dist)./dist;

I_wfs_perfect = zeros(numDist_ps, numDist_rec, numX_rec);
I_wfs = zeros(numDist_ps, numDist_rec, numX_rec);
I_ps = zeros(numDist_ps, numDist_rec, numX_rec);
for ps = 1:numDist_ps
    fprintf('%d/%d\n', ps, numDist_ps);
    % Calculate the field on every receiving point    
    for yrec = 1:numDist_rec
        for xrec = 1:numX_rec
            g = @(x_s) sqrt(r_rec(d_rec(yrec), x_rec(xrec), x_s)./(r_rec(d_rec(yrec), x_rec(xrec), x_s) + r_ps(d_ps(ps), x_s)));
            G_ = @(x_s) G(r_rec(d_rec(yrec), x_rec(xrec), x_s), k);
            
            Q_perfect_ = @(x_s) Q_perfect(x_s, d_ps(ps));
            aux = @(x_s) Q_perfect_(x_s).*G_(x_s).*g(x_s);
            I_wfs_perfect(ps, yrec, xrec) = integral(aux, limx(1), limx(2));
            
            Q_ = @(x_s) Q(x_s, d_ps(ps));
            aux = @(x_s) Q_(x_s).*G_(x_s).*g(x_s);
            I_wfs(ps, yrec, xrec) = integral(aux, limx(1), limx(2));
            
            dist = sqrt((d_rec(yrec) + d_ps(ps))^2 + x_rec(xrec)^2);
            I_ps(ps, yrec, xrec) = exp(-1i*k*dist)/dist;
        end
    end
end

repAbs = cat(4, 20*log10(abs(I_wfs_perfect)./abs(I_ps)), 20*log10(abs(I_wfs)./abs(I_ps)));
repPha = cat(4, rad2deg(angle(I_wfs_perfect./I_ps)), rad2deg(angle(I_wfs./I_ps)));

% visualObj = animation({d_ps, d_rec, x_rec, 1:2},...
%     {repAbs}, {'Distance primary source', 'Distance receiver', 'Horiz. deviation receiver', 'Type'}, ...
%     {'|Prel| (dB)'}, [], []);
visualObj = animation({d_ps, d_rec, x_rec, 1:2},...
    {repPha}, {'Distance primary source', 'Distance receiver', 'Horiz. deviation receiver', 'Type'}, ...
    {'|Prel| (dB)'}, [], []);

vecs = {d_ps, d_rec, x_rec, [1 2]};

% Primary source and receiver near field
ax = drawArray(vecs, repAbs, [1, 2], 'indepDimLimits', [0 3; 0 3], 'nonIndepDimIndices', [1 2], 'surf', true);
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.ZLabel.Interpreter = 'latex';
ax.XLabel.String = '$d_{ps}$';
ax.YLabel.String = '$d$';
ax.ZLabel.String = '$20\log_{10}(|P_{rel}|)$';

% Receiver far field, primary source near field
ax = drawArray(vecs, repAbs, 1, 'indepDimLimits', [0 3], 'nonIndepDimValues', [10 0 2]);
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.XLabel.String = '$d_{ps}/\lambda$';
ax.YLabel.String = '$20\log_{10}(|P_{rel}|)$ (dB)';
ax.XLabel.FontSize = 15;
ax.YLabel.FontSize = 15;
[aux, x, y] = drawArray(vecs, repPha, 1, 'indepDimLimits', [0 3], 'nonIndepDimValues', [Inf 0 2]);
close(aux.Parent)
yyaxis(ax, 'right')
plot(ax, x, y)
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = '$\angle{P_{rel}} (^{\circ})$';
ax.YLabel.FontSize = 15;
% printfig(ax.Parent, imagesPath, 'Experiment10_TruncationIdealCaseA', 'eps')

% Primary source far field, receiver near field
ax = drawArray(vecs, repAbs, 2, 'indepDimLimits', [0 3], 'nonIndepDimValues', [Inf 0 2]);
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.XLabel.String = '$d/\lambda$';
ax.YLabel.String = '$20\log_{10}(|P_{rel}|)$ (dB)';
ax.XLabel.FontSize = 15;
ax.YLabel.FontSize = 15;
[aux, x, y] = drawArray(vecs, repPha, 2, 'indepDimLimits', [0 3], 'nonIndepDimValues', [Inf 0 2]);
close(aux.Parent)
yyaxis(ax, 'right')
plot(ax, x, y)
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = '$\angle{P_{rel}} (^{\circ})$';
ax.YLabel.FontSize = 15;
% printfig(ax.Parent, imagesPath, 'Experiment10_TruncationIdealCaseB', 'eps')

% Observation: I don't see a significative difference in performance. So
% don't worry too much about the 1/r0 term, because the dominant limitation
% is the inaccuracy of the stationary phase method

% B) d_ps = inf
Lmax = 100;
numDist_rec = 10;
d_rec = linspace(1e-1, 100, numDist_rec);
d_ps = 100;
k = 1;

dl = 0.01;
x_s = 0:dl:Lmax/2;
numL = length(x_s);

I_wfs = zeros(numDist_rec, numL);
for yrec = 1:numDist_rec
    g = @(x_s) sqrt(r_rec(d_rec(yrec), 0, x_s)./(r_rec(d_rec(yrec), 0, x_s) + r_ps(d_ps, x_s)));
    G_ = @(x_s) G(r_rec(d_rec(yrec), 0, x_s), k);
    
    Q_ = @(x_s) Q(x_s, d_ps);
    aux = @(x_s) Q_(x_s).*G_(x_s).*g(x_s);
    
    val = aux(x_s);
    val(2:end) = val(2:end)*2;
    I_wfs(yrec, :) = cumsum(val)*dl;
end

I_ps = (exp(-1i*k*(d_ps + d_rec))./(d_ps + d_rec)).';

I_psExt = pointWiseExtend(I_ps, I_wfs);

repAbs = 20*log10(abs(I_wfs(:, 1:50:end))./abs(I_psExt(:, 1:50:end)));
visualObj = animation({d_rec, x_s(1:50:end)*2},...
    {repAbs}, {'Distance receiver', 'L (m)'}, ...
    {'|Prel| (dB)'}, [], []);

% Other way more convenient for plotting
L = 1;
numDist_rec = 100;
numDist_ps = 1;
d_rec = linspace(0, 10, numDist_rec+1); d_rec = d_rec(2:end);
d_ps = 100;
k = 2*pi*(0.5:0.2:50);
numK = length(k);

Q = @(x_s, d_ps, k) d_ps./r_ps(d_ps, x_s).*exp(-1i*k*r_ps(d_ps, x_s))./sqrt(r_ps(d_ps, x_s)).*sqrt(1i*k/(2*pi)); % Feeding of each infinitesimal secondary source

I_wfs = zeros(numDist_rec, numK);
I_ps = zeros(numDist_rec, numK);
for yrec = 1:numDist_rec
    fprintf('%d/%d\n', yrec, numDist_rec)
    g = @(x_s) sqrt(r_rec(d_rec(yrec), 0, x_s)./(r_rec(d_rec(yrec), 0, x_s) + r_ps(d_ps, x_s)));
    for kIter = 1:numK
        G_ = @(x_s) G(r_rec(d_rec(yrec), 0, x_s), k(kIter));
        Q_ = @(x_s) Q(x_s, d_ps, k(kIter));
        aux = @(x_s) Q_(x_s).*G_(x_s).*g(x_s);
        
        dl = min(0.01/k(kIter), 0.01);
        x_s = 0:dl:L/2;
        
        val = aux(x_s);
        val(2:end) = val(2:end)*2;
                
        I_wfs(yrec, kIter) = sum(val)*dl;
        I_ps(yrec, kIter) = G(d_rec(yrec) + d_ps, k(kIter));
    end
end

repAbs = 20*log10(abs(I_wfs)./abs(I_ps));
repPha = rad2deg(angle(I_wfs./I_ps));

visualObj = animation({d_rec, k},...
    {repAbs}, {'Distance receiver', 'k (rad/m)'}, ...
    {'|Prel| (dB)'}, [], []);

visualObj = animation({d_rec, k},...
    {repPha}, {'Distance receiver', 'k (rad/m)'}, ...
    {'\angle{Prel} (\circ)'}, [], []);

% Find the k for which the angle starts to converge
[minVal, minInd] = min(repPha, [], 2);
kThreshold = k(minInd);
ax = axes(figure);
plot(ax, d_rec, kThreshold)

ax = drawArray(vecs, repAbs, [2 1], 'indepDimLimits', [0 Inf], 'nonIndepDimValues', [Inf 0 2]);

% C) Line source integral
f = @(x_s, d, k) exp(-1i*k*sqrt(d^2 + x_s.^2))./sqrt(d^2 + x_s.^2);

d = 10;
k = 1;
Lmax = 100;
dx = 2*pi/k*0.01;
sampVec = 0:dx:Lmax/2;
val = f(sampVec, d, k)*dx;
val(2:end) = 2*val(2:end);

acumI = cumsum(val);
ideal = exp(-1i*k*d)*sqrt(2*pi/(1i*k*d)); % When we have integrated to the infinity

cmap = lines;
ax = axes(figure, 'NextPlot', 'Add');
plot(ax, sampVec*2, abs(acumI), 'Color', cmap(1,:));
plot(ax, [0, sampVec(end)*2], abs(ideal)*ones(1, 2), '--', 'Color', cmap(1,:))
ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = 'D (m)';
ax.YLabel.String = '$|I|$';
ax.YLabel.Interpreter = 'latex';
ax.YColor = cmap(1,:);
yyaxis(ax, 'right')
plot(ax, sampVec*2, rad2deg(angle(acumI)), 'Color', cmap(2,:));
plot(ax, [0, sampVec(end)*2], rad2deg(angle(ideal))*ones(1, 2), '--', 'Color', cmap(2,:))
ax.YLabel.String = '$\angle{I} (^\circ)$';
ax.YLabel.Interpreter = 'latex';
ax.YColor = cmap(2,:);
yyaxis('left')
ax.YColor = cmap(1,:);
ax.Title.String = '$d = 10$, $k = 1$';
ax.Title.Interpreter = 'latex';

% printfig(ax.Parent, imagesPath, 'lineSourceIntegral', 'eps');

% Not very useful. It's better to normalize:
fnorm = @(x_s, lambdaNorm, dNorm) exp(-1i*2*pi/lambdaNorm*sqrt(dNorm^2 + x_s.^2))./sqrt(dNorm^2 + x_s.^2);
% lambdaNorm = lambda/L and dNorm = d/L

dNorm = linspace(0, 2, 100); dNorm = dNorm(2:end);
LNorm = linspace(0, 20, 1000); LNorm = LNorm(2:end);

numD = length(dNorm);
numL = length(LNorm);

I = zeros(numD, numL);
for dIter = 1:numD
    for Liter = 1:numL
        lambdaNorm = 1/LNorm(Liter);
        dx = min(lambdaNorm*0.01, 0.01);
        sampVec = 0:dx:1/2;
        val = fnorm(sampVec, lambdaNorm, dNorm(dIter))*dx;
        val(2:end) = 2*val(2:end);
        
        ideal = exp(-1i*2*pi/lambdaNorm*dNorm(dIter))*sqrt(2*pi/(1i*2*pi/lambdaNorm*dNorm(dIter)));
        I(dIter, Liter) = sum(val)/ideal;
    end
end

visualObj = animation({dNorm, LNorm},...
    {rad2deg(angle(I))}, {'d/L', 'L/\lambda'}, ...
    {'angle(I)'}, [], []);

% Find the k for which the angle starts to converge
[minVal, minInd] = min(rad2deg(angle(I)), [], 2);
LNormThreshold = LNorm(minInd);
% When d/L (dNorm) is very small, the method of finding the minimum is not useful
% because the shape of the curve is different. Let's remove those values
% when representing. To do that, we put the condition that LNormThreshold
% values should be monotically increasing with dNorm
validInd = (LNormThreshold(2:end) - LNormThreshold(1:end-1)) >= 0;
ini = find(validInd, 1, 'first');
ax = axes(figure);
plot(ax, dNorm(ini:end), LNormThreshold(ini:end))
ax.XLabel.String = '$d/D$';
ax.YLabel.String = '$D/\lambda$';
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.YLim = [0, 15];
% printfig(ax.Parent, imagesPath, 'LnormThreshold', 'eps');

% GTAC side
s = WFSToolSimple.generateScenario(96);
Lx = max(s.loudspeakersPosition(:,1)) - min(s.loudspeakersPosition(:,1));
L = 0.18*23;
c = 340;
dNormGTAC = linspace(0, Lx, 100)/L;
f = 1:0.25:500;

numD = length(dNormGTAC);
numF = length(f);
I = zeros(numD, numF);
for dIter = 1:numD
    for fIter = 1:numF
        lambdaNorm = c/(f(fIter)*L);
        dx = min(lambdaNorm*0.01, 0.01);
        sampVec = 0:dx:1/2;
        val = fnorm(sampVec, lambdaNorm, dNormGTAC(dIter))*dx;
        val(2:end) = 2*val(2:end);
        
        ideal = exp(-1i*2*pi/lambdaNorm*dNormGTAC(dIter))*sqrt(2*pi/(1i*2*pi/lambdaNorm*dNormGTAC(dIter)));
        I(dIter, fIter) = sum(val)/ideal;
    end
end

visualObj = animation({dNormGTAC, f},...
    {rad2deg(angle(I))}, {'d/L', 'f'}, ...
    {'angle(I)'}, [], []);

% Find the k for which the angle starts to converge
[minVal, minInd] = min(rad2deg(angle(I)), [], 2);
fThreshold = f(minInd);
% When d/L (dNorm) is very small, the method of finding the minimum is not useful
% because the shape of the curve is different. Let's remove those values
% when representing. To do that, we put the condition that LNormThreshold
% values should be monotically increasing with dNorm
notValidInd = (fThreshold(2:end) - fThreshold(1:end-1)) < 0;
ini = find(notValidInd, 1, 'last') + 1;
ax = axes(figure);
plot(ax, dNormGTAC(ini:end)*L, fThreshold(ini:end))
ax.XLabel.String = '$d$ (m)';
ax.YLabel.String = '$f$ (Hz)';
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
% printfig(ax.Parent, imagesPath, 'freqThresholdGTAC', 'eps');

%% Scheme of the truncation scenario

ax = axes(figure, 'NextPlot', 'Add');
D = 1;
d_ps = 1;
d = 0.5;
plot(ax, [-D/2, D/2], [0 0], 'k');
scatter(ax, 0, -d_ps, 50, [0 0 0], 'filled');
a = scatter(ax, 0, d, 50, [0 0 0], 'filled');

drawbrace([0 0], [0 -d_ps], 10); % Add curly brace
drawbrace([0 0], [0 d], 10); % Add curly brace
drawbrace([0 0], [D/2 0], 10); % Add curly brace
text(ax, 0.05, -d_ps/2, '$d_{\PosTheoSubInd[primarySource]}$');
text(ax, -0.1, d/2, '$d$');
text(ax, D/4, 0.1, '$\sectionTheoLength/2$');
text(ax, 0.05, -d_ps, '$\PosTheo[primarySource]$')
text(ax, 0.05, d, '$\PosTheo$')
% figCoord = data2figureCoordinates( ax, [-D/2*0.9, -d_ps*0.1; -D/2*0.9, 0] );
% annotation('textarrow', [figCoord(1, 1), figCoord(2, 1)], [figCoord(1, 2), figCoord(2, 2)], 'String', '$\sectionTheo$')
text(ax, -D/2*0.9, -d_ps*0.1, '$\sectionTheo$')

ax.Visible = 'off';

% currentFolder = pwd;
% cd(imagesPath); % Needed for inkscape to link svg files properly
% Plot2LaTeX( ax.Parent, 'Experiment10_truncationScheme'); % 'Experiment10_truncationScheme'
% cd(currentFolder)
