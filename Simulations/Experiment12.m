%% Experiment 12

% Created on 2nd July 2018.
% It generates the graphs for the section Compromise between frequency
% filter length and accuracy.

%% Preamble
pathSetUp;

imagesPath = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\TFM\Img\';

dataPathName = [globalPath, 'Data\'];
ID = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

%% How does sampling frequency, number of coefficients and performance of frequency filter relate?

fs = 44100;
c = 340;

% Filter variables for the time WFS filter. Creation of frequency filters
% with different orders.
magnFiltOrder = 2.^(8);
hilbertFiltOrder = 2.^(10);
fs = 44100*[0.5];
freqFilters = cell(numel(fs), 1);
hMags = cell(numel(fs), 1);
hHilbs = cell(numel(fs), 1);
for k = 1:numel(fs)
[freqFilter, freqFiltDelay, hMag, delayMag, hHilb, delayHilb] = getFrequencyFilter( magnFiltOrder, hilbertFiltOrder, fs(k) );
freqFilters{k} = freqFilter;
hMags{k} = hMag;
hHilbs{k} = hHilb;
end

ax = axes(figure, 'NextPlot', 'Add');
ax.Title.String = 'H(f)';
ax.XLim = [-1000, 1000];
ax.XLabel.String = 'Frequency (Hz)';
ax.YLim = [0, 4];
ax.YLabel.String = '|H(f)|';
yyaxis(ax, 'right')
ax.YLim = [-100, 100];
ax.YLabel.String = '$\angle{H(f)}$'; ax.YLabel.Interpreter = 'latex';

axMag = axes(figure, 'NextPlot', 'Add');
axMag.Title.String = 'H_{mag}';
axMag.XLim = [-1000, 1000];
axMag.XLabel.String = 'Frequency (Hz)';
axMag.YLim = [0, 4];
axMag.YLabel.String = '|H_{mag}|';
yyaxis(axMag, 'right')
axMag.YLim = [-180 180];
axMag.YLabel.String = '$\angle{H_{mag}}$'; axMag.YLabel.Interpreter = 'latex';

axHilb = axes(figure, 'NextPlot', 'Add');
axHilb.XLim = [-1000, 1000];
axHilb.XLabel.String = 'Frequency (Hz)';
axHilb.YLim = [0 2];
axHilb.YLabel.String = '|H_{Hilb}(f)|';
yyaxis(axHilb, 'right');
axHilb.YLim = [-100 100];
axHilb.YLabel.String = '$\angle{H_{Hilb}(f)}$'; axHilb.YLabel.Interpreter = 'latex';
axHilb.Title.String = 'H_{Hilb}';

for k = 1:numel(fs)
    f = -1000:1000;
    
    % Total magnitude
    freqFilter = freqFilters{k};
    Htotal = freqz(freqFilter, 1, f, fs(k));
    yyaxis(ax, 'left')
    plot(ax, f, abs(Htotal))
    yyaxis(ax, 'right')
    phaseShift = -freqFiltDelay/fs(k)*2*pi*f;
    pha = rad2deg(angle(Htotal.*exp(-1i*phaseShift)));
    plot(ax, f, pha) 
           
    % Magnitude Filter
    hMag = hMags{k};
    HMag = freqz(hMag, 1, f, fs(k));
    yyaxis(axMag, 'left');
    plot(axMag, f, abs(HMag))
    yyaxis(axMag, 'right');
    phaseShift = -delayMag/fs(k)*2*pi*f;
    plot(axMag, f, rad2deg(angle(HMag.*exp(-1i*phaseShift))))
    
    % Hilbert
    hHilb = hHilbs{k};
    HHilb = freqz(hHilb, 1, f, fs(k));
    yyaxis(axHilb, 'left');
    plot(axHilb, f, abs(HHilb))
    yyaxis(axHilb, 'right');
    phaseShift = -delayHilb/fs(k)*2*pi*f;
    plot(axHilb, f, rad2deg(angle(HHilb.*exp(-1i*phaseShift))))
end
yyaxis(ax, 'left')
plot(ax, f, sqrt(abs(f)/c))
yyaxis(axMag, 'left')
plot(axMag, f, sqrt(abs(f)/c))

str = cell(numel(fs), 1);
for k = 1:numel(fs)
    str{k} = ['fs = ', num2str(fs(k))];
end
legend(ax, [str; {'Ideal'}])
legend(axHilb, str)
legend(axMag, [str; {'Ideal'}])

% 
% ax = axes(figure, 'NextPlot', 'Add');
% for k = 1:numel(fs)
% t = ((0:N-1) - freqFiltDelay)/fs(k);
% plot(ax, t, freqFilters{k})
% end

% Conclussion: what matters in terms of performance is the duration of the
% filter in seconds. As long as the sampling frequency is bigger than the
% nyquist frequency (2*fmax), the duration of the FIR filter in seconds
% (numberOfCoefficients/fs) is what determines the accuracy.
% Now, the question is: how long should it be in order to have an
% acceptable accuracy? We will test both filters separately.

%% Test magnitude filter. How does number of coefficients and performance relate?
magnFiltOrder = 2.^(2:9);
Nmag = length(magnFiltOrder);
hilbertFiltOrder = 2^14;
fs = 44100;
c = 340;
hMags = cell(Nmag, 1);
delaysMag = zeros(Nmag, 1);
for k = 1:Nmag
    [hMag, delayMag] = getFrequencyFilter( magnFiltOrder(k), [], fs );
    hMags{k} = hMag;
    delaysMag(k) = delayMag;
end
[hHilbert, delayHilbert] = getFrequencyFilter( [], hilbertFiltOrder, fs );
delta = zeros(size(hHilbert));
delta(delayHilbert + 1) = 1;

hTotals = cell(Nmag, 1);
delays = zeros(Nmag, 1);
for n = 1:Nmag  
    hTotal = 1/sqrt(2)*conv((delta - hHilbert), hMags{n});
    delay = delaysMag(n) + delayHilbert;
    
    hTotals{n} = hTotal;
    delays(n) = delay;
end

% Filter variables for the time WFS filter. Creation of frequency filters
% with different orders.  
freqFilters = hTotals;
freqFiltDelays = delays;

fRep = -1000:1000;

axHilb = axes(figure, 'NextPlot', 'Add');
axHilb.XLim = [-1000, 1000];
axHilb.XLabel.String = 'Frequency (Hz)';
axHilb.YLim = [0 2];
axHilb.YLabel.String = '|H_{Hilb}(f)|';
yyaxis(axHilb, 'right');
axHilb.YLim = [-100 100];
axHilb.YLabel.String = '$\angle{H_{Hilb}(f)}$'; axHilb.YLabel.Interpreter = 'latex';
axHilb.Title.String = 'H_{Hilb}';
% Hilbert
HHilb = freqz(hHilbert, 1, fRep, fs);
yyaxis(axHilb, 'left');
plot(axHilb, fRep, abs(HHilb))
yyaxis(axHilb, 'right');
phaseShift = -delayHilbert/fs*2*pi*fRep;
plot(axHilb, fRep, rad2deg(angle(HHilb.*exp(-1i*phaseShift))))

ax = axes(figure, 'NextPlot', 'Add');
ax.Title.String = 'H(f)';
ax.XLim = [-1000, 1000];
ax.XLabel.String = 'Frequency (Hz)';
ax.YLim = [0, 4];
ax.YLabel.String = '|H(f)|';
yyaxis(ax, 'right')
ax.YLim = [-100, 100];
ax.YLabel.String = '$\angle{H(f)}$'; ax.YLabel.Interpreter = 'latex';

axMag = axes(figure, 'NextPlot', 'Add');
axMag.Title.String = 'H_{mag}';
axMag.XLim = [-1000, 1000];
axMag.XLabel.String = 'Frequency (Hz)';
axMag.YLim = [0, 4];
axMag.YLabel.String = '|H_{mag}|';
yyaxis(axMag, 'right')
axMag.YLim = [-180 180];
axMag.YLabel.String = '$\angle{H_{mag}}$'; axMag.YLabel.Interpreter = 'latex';

lMagsAbs = gobjects(Nmag, 1);
lMagsPha = gobjects(Nmag, 1);
lTotalAbs = gobjects(Nmag, 1);
lTotalPha = gobjects(Nmag, 1);
for n = 1:Nmag
    % Magnitude Filter
    hMag = hMags{n};
    HMag = freqz(hMag, 1, fRep, fs);
    yyaxis(axMag, 'left');
    lMagsAbs(n) = plot(axMag, fRep, abs(HMag));
    yyaxis(axMag, 'right');
    phaseShift = -delaysMag(n)/fs*2*pi*fRep;
    lMagsPha(n) = plot(axMag, fRep, rad2deg(angle(HMag.*exp(-1i*phaseShift))));
    
    % Total magnitude
    Htotal = freqz(hTotals{n}, 1, fRep, fs);
    yyaxis(ax, 'left')
    lTotalAbs(n) = plot(ax, fRep, abs(Htotal));
    yyaxis(ax, 'right')
    phaseShift = -delays(n)/fs*2*pi*fRep;
    pha = rad2deg(angle(Htotal.*exp(-1i*phaseShift)));
    lTotalPha(n) = plot(ax, fRep, pha);
end
yyaxis(ax, 'left')
plot(ax, fRep, sqrt(abs(fRep)/c));
yyaxis(axMag, 'left')
plot(axMag, fRep, sqrt(abs(fRep)/c))

% aux = repmat({'-'}, Nmag, 1);
% [lMagsAbs(1:Nmag).LineStyle] = aux{:};
cmap = flipud(colorcube);
color1 = cmap(9, :);
color2 = cmap(14, :);
cmapAbs = interp1([1,Nmag], [color1; color2], 1:Nmag);
color1 = cmap(21, :);
color2 = cmap(26, :);
cmapPha = interp1([1,Nmag], [color1; color2], 1:Nmag);
for n = 1:Nmag
    lMagsAbs(n).Color = cmapAbs(n,:);
    lMagsAbs(n).Marker = 'none';
    lMagsAbs(n).LineStyle = '-';
    
    lMagsPha(n).Color = cmapPha(n,:);
    lMagsPha(n).Marker = 'none';
    lMagsPha(n).LineStyle = '-';
    
    lTotalAbs(n).Color = cmapAbs(n,:);
    lTotalAbs(n).Marker = 'none';
    lTotalAbs(n).LineStyle = '-';
    
    lTotalPha(n).Color = cmapPha(n,:);
    lTotalPha(n).Marker = 'none';
    lTotalPha(n).LineStyle = '-';
end

str = cell(Nmag, 1);
for n = 1:Nmag
    str{n} = ['N = ', num2str(magnFiltOrder(n))];
end
legend(ax, [str; {'Ideal'}])
legend(axMag, [str; {'Ideal'}])

%%% Simulation
if ~exist('obj', 'var') || ~isvalid(obj)
    obj = SimulationController;
end

% Constants
WFSfilterLength = 22050;
zPos = 1.65;
WFSarrayOffset = [0.46 2.21 zPos]; % [x, y, z] coordinates. Useful for generating acoustic path IR.
roomDim = [4.48, 9.13, 2.64];

% Noise source coefficient
amplitude = 1;
phase = 0;

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
numPointsPerArc = 4;
radius = [3.6 4 4.4 4.8];
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
NSpositions = [x(:), y(:), zeros(numel(x), 1)];

% Frequencies
freqs = 0:10:1000;

% % Signal
% durSign = 1; % Duration of tone for time processing
% t = (0:ceil(durSign*fs)-1)/fs;
% NSsignal = chirp(t, 20, durSign, 940);
% predefSignals = true;

% Room characteristics and impulse response of chamber
beta = 0; % Average reflection coefficient of the walls of the chamber
WFS_AcPath_previously_calculated = true;
NS_AcPath_previously_calculated = true;
appendFreeSpaceAcPaths = false;

% WFS options
frequencyCorrection = true;
attenuationType = 'Ruben';

% Simulation options
timeDomainActive = true;
fakeTimeProcessing = true;
frequencyDomainActive = false;
automaticLengthModification = false;
saveSignals = true;

SetupParametersScript
AcousticPathCalculationScript
simulationScript

%%% Analysis
% recNS_signalsAux = pointWiseExtend( recNS_signals, rec_signals );
% recWFS_signals = rec_signals - recNS_signalsAux;
% 
% f = 0:1000;
% oper = @(x) freqz(x, 1, f, fs);
% recNS = oneDimOperOverMultiDimArray( oper, recNS_signals, 2);
% recWFS = oneDimOperOverMultiDimArray( oper, recWFS_signals, 2);

% % In order to calculate the global cancellation, it is convenient to create
% % a structure array
% ss = repmat(s(1), [numel(fsel), numNSpos]);
% for fIter = 1:numel(fsel)
%     for ns = 1:numNSpos
%         ss(fIter, ns).recCoef = recNS(:, fIter, 1, ns) + recWFS(:, fIter, 1, ns);
%         ss(fIter, ns).recNScoef = recNS(:, fIter, 1, ns); 
%         ss(fIter, ns).recWFScoef = recWFS(:, fIter, 1, ns);
%         ss(fIter, ns).Frequency = fsel(fIter);
%     end
% end

[sExt, corrFactInd, corrFactGlob, attenInd, attenGlob, corrFactAver, gainAver] =...
    SimulationController.addCancellationParametersToStructure(s);

freqEdges = 0:20:1000;
attenEdges = -20:5;
vecs = {magnFiltOrder, freqs, 1:numNSpos};
axCanc = gobjects(numFreqFilters, 1);
for k = 1:numFreqFilters
    [~, attenGlobCurrent] = filterArrayForRepresentation(vecs, attenGlob, [2 3], 'nonIndepDimIndices', k);
    ax = histogram2D(10*log10(attenGlobCurrent), 1, freqs, freqEdges, attenEdges);
    ax.XLabel.String = 'Frequency (Hz)';
    ax.YLabel.String = 'Attenuation (dB)';
    ax.Title.String = ['N = ', num2str(magnFiltOrder(k))];
    colorbar(ax);
    axCanc(k) = ax;
end

% % Print
% sel = 2:5;
% for k = 1:numel(sel)
%     printfig(axCanc(sel(k)).Parent, imagesPath, ['Experiment12_globalAttenMagnOrder_', num2str(magnFiltOrder(sel(k)))], 'eps')
% end
% delaysMag/fs*c

% Average gain
freqEdges = 0:20:1000;
gainEdges = -20:5;
vecs = {magnFiltOrder, freqs, 1:numNSpos};
axGainAver = gobjects(numFreqFilters, 1);
for k = 1:numFreqFilters
    [~, gainAverCurrent] = filterArrayForRepresentation(vecs, gainAver, [2 3], 'nonIndepDimIndices', k);
    ax = histogram2D(10*log10(gainAverCurrent), 1, freqs, freqEdges, gainEdges);
    ax.XLabel.String = 'Frequency (Hz)';
    ax.YLabel.String = 'Average Gain (dB)';
    ax.Title.String = ['N = ', num2str(magnFiltOrder(k))];
    colorbar(ax);
    axGainAver(k) = ax;
end

% % Print
sel = 2:8;
for k = 1:numel(sel)
    printfig(axGainAver(sel(k)).Parent, imagesPath, ['Experiment12_GainAverMagnOrder_', num2str(magnFiltOrder(sel(k)))], 'eps')
end
delaysMag/fs*c

%% Test Hilbert filter. How does number of coefficients and performance relate?
magnFiltOrder = 2^10;
hilbertFiltOrder = 2.^(2:13);
Nhil = length(hilbertFiltOrder);
fs = 44100;
c = 340;
hHilbs = cell(Nhil, 1);
delaysHilb = zeros(Nhil, 1);
for n = 1:Nhil
    [hHilbert, delayHilbert] = getFrequencyFilter( [], hilbertFiltOrder(n), fs );
    hHilbs{n} = hHilbert;
    delaysHilb(n) = delayHilbert;
end
% [hMag, delayMag] = getFrequencyFilter( magnFiltOrder, [], fs );

hTotals = cell(Nhil, 1);
delays = zeros(Nhil, 1);
for n = 1:Nhil
    hHilbert = hHilbs{n};
    delta = zeros(size(hHilbert));
    delta(delaysHilb(n) + 1) = 1;
    hTotal = 1/sqrt(2)*conv((delta - hHilbert), hMag);
    delay = delaysHilb(n) + delayMag;
    
    hTotals{n} = hTotal;
    delays(n) = delay;
end

% Filter variables for the time WFS filter. Creation of frequency filters
% with different orders.  
freqFilters = hTotals;
freqFiltDelays = delays;

% Visualize filters response
fRep = -1000:1000;

axMag = axes(figure, 'NextPlot', 'Add');
axMag.Title.String = 'H_{mag}';
axMag.XLim = [-1000, 1000];
axMag.XLabel.String = 'Frequency (Hz)';
axMag.YLim = [0, 4];
axMag.YLabel.String = '|H_{mag}|';
yyaxis(axMag, 'right')
axMag.YLim = [-180 180];
axMag.YLabel.String = '$\angle{H_{mag}}$'; axMag.YLabel.Interpreter = 'latex';
% Magnitude Filter
HMag = freqz(hMag, 1, fRep, fs);
yyaxis(axMag, 'left');
plot(axMag, fRep, abs(HMag));
yyaxis(axMag, 'right');
phaseShift = -delayMag/fs*2*pi*fRep;
plot(axMag, fRep, rad2deg(angle(HMag.*exp(-1i*phaseShift))));

axHilb = axes(figure, 'NextPlot', 'Add');
axHilb.XLim = [-1000, 1000];
axHilb.XLabel.String = 'Frequency (Hz)';
axHilb.YLim = [0 2];
axHilb.YLabel.String = '|H_{Hilb}(f)|';
yyaxis(axHilb, 'right');
axHilb.YLim = [-100 100];
axHilb.YLabel.String = '$\angle{H_{Hilb}(f)}$'; axHilb.YLabel.Interpreter = 'latex';
axHilb.Title.String = 'H_{Hilb}';

ax = axes(figure, 'NextPlot', 'Add');
ax.Title.String = 'H(f)';
ax.XLim = [-1000, 1000];
ax.XLabel.String = 'Frequency (Hz)';
ax.YLim = [0, 4];
ax.YLabel.String = '|H(f)|';
yyaxis(ax, 'right')
ax.YLim = [-100, 100];
ax.YLabel.String = '$\angle{H(f)}$'; ax.YLabel.Interpreter = 'latex';

lHilbAbs = gobjects(Nmag, 1);
lHilbPha = gobjects(Nmag, 1);
lTotalAbs = gobjects(Nmag, 1);
lTotalPha = gobjects(Nmag, 1);
for n = 1:Nhil
    % Hilbert
    HHilb = freqz(hHilbs{n}, 1, fRep, fs);
    yyaxis(axHilb, 'left');
    lHilbAbs(n) = plot(axHilb, fRep, abs(HHilb));
    yyaxis(axHilb, 'right');
    phaseShift = -delaysHilb(n)/fs*2*pi*fRep;
    lHilbPha(n) = plot(axHilb, fRep, rad2deg(angle(HHilb.*exp(-1i*phaseShift))));
    
    % Total magnitude
    Htotal = freqz(hTotals{n}, 1, fRep, fs);
    yyaxis(ax, 'left')
    lTotalAbs(n) = plot(ax, fRep, abs(Htotal));
    yyaxis(ax, 'right')
    phaseShift = -delays(n)/fs*2*pi*fRep;
    pha = rad2deg(angle(Htotal.*exp(-1i*phaseShift)));
    lTotalPha(n) = plot(ax, fRep, pha);
end

% aux = repmat({'-'}, Nmag, 1);
% [lMagsAbs(1:Nmag).LineStyle] = aux{:};
cmap = flipud(colorcube);
color1 = cmap(9, :);
color2 = cmap(14, :);
cmapAbs = interp1([1,Nhil], [color1; color2], 1:Nhil);
color1 = cmap(21, :);
color2 = cmap(26, :);
cmapPha = interp1([1,Nhil], [color1; color2], 1:Nhil);
for n = 1:Nhil
    lHilbAbs(n).Color = cmapAbs(n,:);
    lHilbAbs(n).Marker = 'none';
    lHilbAbs(n).LineStyle = '-';
    
    lHilbPha(n).Color = cmapPha(n,:);
    lHilbPha(n).Marker = 'none';
    lHilbPha(n).LineStyle = '-';
    
    lTotalAbs(n).Color = cmapAbs(n,:);
    lTotalAbs(n).Marker = 'none';
    lTotalAbs(n).LineStyle = '-';
    
    lTotalPha(n).Color = cmapPha(n,:);
    lTotalPha(n).Marker = 'none';
    lTotalPha(n).LineStyle = '-';
end

str = cell(Nhil, 1);
for n = 1:Nhil
    str{n} = ['N = ', num2str(hilbertFiltOrder(n))];
end
legend(ax, str)
legend(axHilb, str)

%%% Simulation %%%
if ~exist('obj', 'var') || ~isvalid(obj)
    obj = SimulationController;
    obj.WFSToolObj.fig.HandleVisibility = 'off';
end

% Constants
WFSfilterLength = 22050;
zPos = 1.65;
WFSarrayOffset = [0.46 2.21 zPos]; % [x, y, z] coordinates. Useful for generating acoustic path IR.
roomDim = [4.48, 9.13, 2.64];

% Noise source coefficient
amplitude = 1;
phase = 0;

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
numPointsPerArc = 4;
radius = [3.6 4 4.4 4.8];
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
NSpositions = [x(:), y(:), zeros(numel(x), 1)];

% Frequencies
freqs = 0:10:1000;

% % Signal
% durSign = 1; % Duration of tone for time processing
% t = (0:ceil(durSign*fs)-1)/fs;
% NSsignal = chirp(t, 20, durSign, 940);
% predefSignals = true;
% saveSignals = true;

% Room characteristics and impulse response of chamber
beta = 0; % Average reflection coefficient of the walls of the chamber
WFS_AcPath_previously_calculated = true;
NS_AcPath_previously_calculated = true;
appendFreeSpaceAcPaths = false;

% WFS options
frequencyCorrection = true;
attenuationType = 'Ruben';

% Simulation options
timeDomainActive = true;
fakeTimeProcessing = true;
frequencyDomainActive = false;
automaticLengthModification = false;

SetupParametersScript
AcousticPathCalculationScript
simulationScript

%%% Analysis
[sExt, corrFactInd, corrFactGlob, attenInd, attenGlob, corrFactAver, averGain] =...
    SimulationController.addCancellationParametersToStructure(s);

freqEdges = 0:20:1000;
attenEdges = -20:5;
vecs = {hilbertFiltOrder, freqs, 1:numNSpos};
axCanc = gobjects(numFreqFilters, 1);
for k = 1:numFreqFilters
    [~, attenGlobCurrent] = filterArrayForRepresentation(vecs, attenGlob, [2 3], 'nonIndepDimIndices', k);
    ax = histogram2D(10*log10(attenGlobCurrent), 1, freqs, freqEdges, attenEdges);
    ax.XLabel.String = 'Frequency (Hz)';
    ax.YLabel.String = 'Global gain (dB)';
    ax.Title.String = ['N = ', num2str(hilbertFiltOrder(k))];
    ax.Parent.Name = ['N = ', num2str(hilbertFiltOrder(k))];
    colorbar(ax);
    axCanc(k) = ax;
end

% % Print
% sel = 8:11;
% for k = 1:numel(sel)
%     printfig(axCanc(sel(k)).Parent, imagesPath, ['Experiment12_globalAttenHilbOrder_', num2str(hilbertFiltOrder(sel(k)))], 'eps')
% end
% delaysHilb/fs*c

freqEdges = 0:20:1000;
attenEdges = -20:5;
vecs = {hilbertFiltOrder, freqs, 1:numNSpos};
axGainAver = gobjects(numFreqFilters, 1);
for k = 1:numFreqFilters
    [~, averGainCurrent] = filterArrayForRepresentation(vecs, averGain, [2 3], 'nonIndepDimIndices', k);
    ax = histogram2D(10*log10(averGainCurrent), 1, freqs, freqEdges, attenEdges);
    ax.XLabel.String = 'Frequency (Hz)';
    ax.YLabel.String = 'Average gain (dB)';
    ax.Title.String = ['N = ', num2str(hilbertFiltOrder(k))];
    ax.Parent.Name = ['N = ', num2str(hilbertFiltOrder(k))];
    colorbar(ax);
    axGainAver(k) = ax;
end

% Print
sel = 8:11;
for k = 1:numel(sel)
    printfig(axGainAver(sel(k)).Parent, imagesPath, ['Experiment12_gainAverHilbOrder_', num2str(hilbertFiltOrder(sel(k)))], 'eps')
end
delaysHilb/fs*c
