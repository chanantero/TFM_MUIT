%% Basic simulation
% Reduced version of the simulations

%% System parameters

% Constants
c = 340; % Sound velocity (m/s)
fs = 44100; % Sample frequency (samples/s)

% Noise source coefficient
amplitude = 1;
phase = 0;

% Creation of frequency filter
magnFiltOrder = 2^(12);
hilbertFiltOrder = 2^(12);
[freqFilter, freqFiltDelay] = getFrequencyFilter( magnFiltOrder, hilbertFiltOrder, fs );
freqFilterLength = length(freqFilter);

% Position of loudspeakers of the WFS array
d = 0.18; % Separation between two contiguous loudspeakers. Size of one loudspeaker
nb = 8; % Bottom and upper sides of the octogon (2 sides)
nd = 8; % Diagonal sides of the octogon (4 sides)
nl = 24; % Lateral side of the octogon (2 sides)
betabd = 45; % Deviation angle between bottom/upper and diagonal sides
[ x, y, alfa ] = octogon(d, nb, nd, nl, betabd);
z = zeros(numel(x), 1);
WFSpositions = [x, y, z];
loudspeakersOrientation = [cosd(alfa), sind(alfa), zeros(numel(alfa), 1)];
numWFS = size(WFSpositions, 1); % 96  
centerX = (max(x) + min(x))/2;
centerY = (max(y) + min(y))/2;

% Microphone position
recPosition = [centerX, centerY, 0];

% Positions of the noise source
NSposition = [centerX + 5, centerY, 0];    
% NSposition = [2.3 0.1 1.5];

durSign = 1; % Duration of signal
t = (0:ceil(durSign*fs)-1)/fs;
NSsignal = chirp(t, 20, durSign, 940);

%% Room characteristics and impulse response of chamber
numReverbTime = 2;
beta = linspace(0, 1, numReverbTime); % Average reflection coefficient of the walls of the chamber
Beta = beta(:) * [1 1 1 1 1 1];
WFSarrayOffset = [0.5, 1.5, 1.5]; % [x, y, z] coordinates. Useful for generating acoustic path IR.
r = recPosition + WFSarrayOffset; % Receiver position [x y z] (m)
wfsPos = WFSpositions + repmat(WFSarrayOffset, numWFS, 1);
nsPos = NSposition + WFSarrayOffset;
maxX = max([nsPos(:, 1); wfsPos(:, 1)]);
maxY = max([nsPos(:, 2); wfsPos(:, 2)]);
roomDim = [maxX+1 maxY+1 4]; % Room dimensions [x y z] (m)

% Adjust the number of samples of impulse responses
dist = sqrt(sum((recPosition - NSposition).^2, 2));
numSampIR = 2^ceil(log2(dist/c*fs)); % Number of samples

WFSfilterLength = numSampIR*2;

WFS_IR = zeros(numWFS, numSampIR, numReverbTime);
for k = 1:numWFS
    for rt = 1:numReverbTime
        WFS_IR(k, :, rt) = rir_generator(c, fs, r, wfsPos(k, :), roomDim, Beta(rt, :), numSampIR);
    end
end

NS_IR = zeros(numReverbTime, numSampIR);
for rt = 1:numReverbTime
    NS_IR(rt, :) = rir_generator(c, fs, r, nsPos, roomDim, Beta(rt, :), numSampIR);
end



%% Simulation
preN = freqFiltDelay;
postN = numSampIR - 1 + freqFiltDelay - 1;
x = [zeros(1, preN), NSsignal, zeros(1, postN)];
numSamp = length(x);

% Generate signals of WFS array loudspeakers

% Calculate delay and attenuation
relPos = WFSpositions - repmat(NSposition, numWFS, 1);
distances = sqrt(sum(relPos.^2, 2));
delays = distances/c;
r0 = 1.44*(0.5+cosd(45));
A = sqrt(r0./(r0 + distances));
cosAlfa = dot(relPos, loudspeakersOrientation, 2)./distances;
attenuations = -A.*cosAlfa./sqrt(distances)*d;
attenuations(cosAlfa < 0) = 0;

% Create basic filter impulse response (delta)
indDelta = floor(delays*fs) + 1;
NsFilt = max(floor(delays(:)*fs)) + 1 + freqFilterLength;
    
filtersWFS = zeros(numWFS, NsFilt);
for ss = 1:numWFS
    if attenuations(ss) ~= 0
        filtersWFS(ss, indDelta(ss)) = attenuations(ss);
    end
end

% Convolute it with the frequency filter
filtersWFS = fftfilt(freqFilter, filtersWFS')';

% Apply it to noise source signal
wfsSignals = fftfilt(filtersWFS', x')';
% Compensate delay
wfsSignals = [wfsSignals(:, freqFiltDelay + 1:end), zeros(numWFS, freqFiltDelay)];

% Simulate
acousticPaths = [WFS_IR; permute(NS_IR, [3 2 1])];
sourceSignals = [wfsSignals; x];

recSignals = zeros(numReverbTime, numSamp);
for rt = 1:numReverbTime
    recSignals(rt, :) = sum(fftfilt(acousticPaths(:, :, rt)', sourceSignals'), 2)';
end

% Simulate with only the noise source
recSignalsNS = zeros(numReverbTime, numSamp);
for rt = 1:numReverbTime
    recSignalsNS(rt, :) = sum(fftfilt(acousticPaths(end, :, rt)', sourceSignals(end,:)'), 2)';
end

%% Visualization
t = (0:numSamp-1)/fs;

axs = gobjects(numReverbTime, 1),
for rt = 1:numReverbTime
    ax = axes(figure);
    plot(ax, t, recSignalsNS(rt,:)', t, recSignals(rt,:)')
    axs(rt) = ax;
end

