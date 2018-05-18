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

% Filter variables for the time WFS filter. 
load('WFSTool/WFSfilter.mat')
freqFilter = hTotal;

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
frequencyCorrection = true;

% Simulation options
timeDomainActive = false;
frequencyDomainActive = true;

%% Setup first parameters

% Noise source variables
obj.amplitude = amplitude;
obj.amplitude(2) = -amplitude;
obj.phase = phase;

% Default values. They don't matter, but do not touch just in case.
obj.NSposition = [3.35 -0.2 0]; % Assumed real position
obj.frequency = 800;
obj.Fs = fs;

% Frequency filter
obj.WFSToolObj.freqFilter = hTotal;
obj.WFSToolObj.filterWFS_length = WFSfilterLength;

% Microphone positions
obj.microPos = recPositions;
numMicro = size(recPositions, 1);

% NS positions
numNSpos = size(NSpositions, 1);

% Room characteristics and impulse responses
Beta = beta(:) * [1 1 1 1 1 1];
numReverbTime = length(beta);

% WFS options
obj.WFSToolObj.frequencyCorrection = frequencyCorrection; % Very important! We want to see what happens without correction

%% Pre-calculate impulse responses
% Since it would take a lot of time to calculate the IR during the
% simulations, we calculate it previously and save it to a .mat. 

r = recPositions + repmat(WFSarrayOffset, numMicro, 1); % Receiver position [x y z] (m)
wfsPos = obj.WFSposition + repmat(WFSarrayOffset, obj.numWFS, 1);
nsPos = NSpositions + repmat(WFSarrayOffset, numNSpos, 1);
maxX = max([nsPos(:, 1); wfsPos(:, 1)]);
maxY = max([nsPos(:, 2); wfsPos(:, 2)]);
roomDim = [maxX+1 maxY+1 4];                % Room dimensions [x y z] (m)

% Minimize the number of samples
dist = zeros(numMicro, numNSpos);
for ns = 1:numNSpos
    dist(:, ns) = sqrt(sum((recPositions - repmat(NSpositions(ns, :), [numMicro, 1])).^2, 2));
end
maxDist = max(dist(:));
numSampIR = 2^ceil(log2(maxDist/c*fs)); % Number of samples

WFSfilterLength = numSampIR*2;
obj.WFSToolObj.filterWFS_length = WFSfilterLength;

if WFS_AcPath_previously_calculated
    WFS_IR = zeros(obj.numMicro, numSampIR, obj.numWFS, numReverbTime);
    disp('IR WFS')
    for k = 1:obj.numWFS
        fprintf('%d/%d\n', k, obj.numWFS);
        disp(' Reverberation Time:')
        for rt = 1:numReverbTime
            fprintf(' %d/%d\n', rt, numReverbTime);
            WFS_IR(:, :, k, rt) = rir_generator(c, fs, r, wfsPos(k, :), roomDim, Beta(rt, :), numSampIR);
        end
    end
    
    % Calculate the frequency responses
    WFS_FR = zeros(obj.numMicro, numFreqs, obj.numWFS, numReverbTime);
    disp('FR WFS')
    for k = 1:obj.numWFS
        fprintf('%d/%d\n', k, obj.numWFS);
        disp(' Reverberation Time:')
        for rt = 1:numReverbTime
            WFS_FR(:, :, k, rt) = DFT_slow(fs*WFS_IR(:, :, k, rt).', fs, freqs).';
        end
    end
    
end

% Do the same with the noise source positions
if NS_AcPath_previously_calculated
    NS_IR = zeros(obj.numMicro, numSampIR, numNSpos, numReverbTime);
    disp('IR NS')
    for k = 1:numNSpos
        fprintf('%d/%d\n', k, numNSpos);
        disp(' Reverberation Time:')
        for rt = 1:numReverbTime
            fprintf(' %d/%d\n', rt, numReverbTime);
            NS_IR(:, :, k, rt) = rir_generator(c, fs, r, nsPos(k,:), roomDim, Beta(rt, :), numSampIR);
        end
    end
    
    NS_FR = zeros(obj.numMicro, numFreqs, numNSpos, numReverbTime);
    disp('FR NS')
    for k = 1:numNSpos
        fprintf('%d/%d\n', k, numNSpos);
        disp(' Reverberation Time:')
        for rt = 1:numReverbTime
            NS_FR(:, :, k, rt) = DFT_slow(fs*NS_IR(:, :, k, rt).', fs, freqs).';
        end
    end
end

if WFS_AcPath_previously_calculated && NS_AcPath_previously_calculated ...
        && appendFreeSpaceAcPaths
    acPath = simulator.calculateMonopolesIR(obj.WFSposition, recPositions, c, fs, numSampIR);
    WFS_IR = cat(4, permute(acPath, [1 3 2 4]), WFS_IR);
    
    acPath = simulator.calculateTheoricAcousticPaths(...
        obj.WFSToolObj.WFSarrayPosition, obj.WFSToolObj.WFSarrayRadiationPattern, obj.WFSToolObj.WFSarrayOrientation,...
        recPositions, obj.WFSToolObj.receiverRadiationPattern, obj.WFSToolObj.receiverOrientation, freqs, c);
    WFS_FR = cat(4, permute(acPath, [1 3 2 4]), WFS_FR);
    
    acPath = simulator.calculateMonopolesIR(NSpositions, recPositions, c, fs, numSampIR);
    NS_IR = cat(4, permute(acPath, [1 3 2 4]), NS_IR);
    
    acPath = simulator.calculateTheoricAcousticPaths(...
        NSpositions, repmat({@simulator.monopoleRadPat}, numNSpos, 1), repmat([1 0 0 1], numNSpos, 1),... % [1 0 0 1] is just a default orientation. It doesn't matter because they are monopoles
        recPositions, obj.WFSToolObj.receiverRadiationPattern, obj.WFSToolObj.receiverOrientation, freqs, c);
    NS_FR = cat(4, permute(acPath, [1 3 2 4]), NS_FR);
    
    numReverbTime = numReverbTime + 1;
end

% WFSposition = obj.WFSposition;
% save([dataPathName, 'AcPathResponsesWFS_', ID, '.mat'], 'r', 'roomDim', 'numSamp', 'WFSposition', 'freqs', 'WFS_IR', 'WFS_FR')
% save([dataPathName, 'AcPathResponsesNS_', ID, '.mat'], 'r', 'roomDim', 'numSamp', 'NSpositions', 'NS_IR', 'freqs', 'NS_FR')

%% Pre-calculate WFS filters

filtersWFS_IR = zeros(obj.numWFS, WFSfilterLength, numNSpos);

obj.domain = 'time';
disp('NS positions:')
for ns = 1:numNSpos
    fprintf('%d/%d\n', ns, numNSpos);
    
    obj.NSposition = NSpositions(ns, :);
    
    obj.WFSToolObj.updateFiltersWFS(); % The NS position has changed, so we need to calculate again the WFS filters
    
    filtersWFS_IR(:, :, ns) = permute(obj.WFSToolObj.filtersWFS_IR(:, 1, :), [1, 3, 2]);
end

%% Simulation
durSign = 1; % Duration of signal
t = (0:ceil(durSign*obj.Fs)-1)/obj.Fs;
[~, filterDelay] = max(obj.WFSToolObj.freqFilter);
filterDelay = filterDelay - 1;

recNScoef_time = zeros(numMicro, numFreqs, numNSpos, numReverbTime);
recCoef_time = zeros(numMicro, numFreqs, numNSpos, numReverbTime);
recWFScoef_time = zeros(numMicro, numFreqs, numNSpos, numReverbTime);
WFScoef_time = zeros(obj.numWFS, numFreqs, numNSpos, numReverbTime);

recNScoef_freq = zeros(obj.numMicro, numFreqs, numNSpos, numReverbTime);
recCoef_freq = zeros(obj.numMicro, numFreqs, numNSpos, numReverbTime);
recWFScoef_freq = zeros(obj.numMicro, numFreqs, numNSpos, numReverbTime);
WFScoef_freq = zeros(obj.numWFS, numFreqs, numNSpos, numReverbTime);

% Progress bar variables
maxIterPerLevel = [numReverbTime, numFreqs, numNSpos];
numLevels = length(maxIterPerLevel);
totalIter = prod(maxIterPerLevel);
levelWeight = flip([1, cumprod(flip(maxIterPerLevel(2:end)))])';
% completedIterPerLevel = [0 0 0]; % Completed
% completedIter = completedIterPerLevel*levelWeight;
% progress = completedIter/totalIter;
wb = waitbar(0, 'Progress');

for rt = 1:numReverbTime
    
    % Set up acoustic paths for the WFS array
    if WFS_AcPath_previously_calculated
        % In case you calculate WFS_IR and WFS_FR previously
        WFSacPathIR = permute(WFS_IR(:, :, :, rt), [1, 3, 2]);
        WFSacPathFR = permute(WFS_FR(:, :, :, rt), [1, 3, 2]);
    else
        WFS_IR_current = zeros(obj.numMicro, numSampIR, obj.numWFS);
        for k = 1:obj.numWFS
            WFS_IR_current(:, :, k) = rir_generator(c, fs, r, wfsPos(k, :), roomDim, Beta(rt, :), numSampIR);
        end
        WFSacPathIR = permute(WFS_IR_current, [1, 3, 2]);
        
        WFS_FR_current = zeros(obj.numMicro, numFreqs, obj.numWFS);
        disp('FR WFS')
        for k = 1:obj.numWFS
            WFS_FR_current(:, :, k) = DFT_slow(fs*WFS_IR_current(:, :, k).', fs, freqs).';
        end
        WFSacPathFR = permute(WFS_FR_current, [1, 3, 2]);
    end
    
    obj.domain= 'time';
    obj.setAcousticPaths('WFS', WFSacPathIR);
    
    obj.domain= 'frequency';
    WFSacPathFRstruct = struct('acousticPaths', WFSacPathFR, 'frequencies', freqs);
    obj.setAcousticPaths('WFS', WFSacPathFRstruct); % obj.setAcousticPaths('WFS', 'theoretical');
    
    for f = 1:numFreqs
        fcurr = freqs(f); % Current frequency
        obj.frequency = fcurr;
        
        % Generate tone with the propper frequency
        x = obj.amplitude(1) * cos(2*pi*fcurr*t + obj.phase(1));
        x = [zeros(1, filterDelay), x, zeros(1, numSampIR - 1)];
               
        obj.domain = 'time';
        obj.NScoef = x;
        obj.NSVcoef = -x;
        
        for ns = 1:numNSpos
            obj.NSposition = NSpositions(ns, :);
            
            % Set up acoustic paths for the noise source
            obj.domain= 'time';
            NSacPathIR = repmat(permute(NS_IR(:, :, ns, rt), [1, 3, 2]), [1, 2, 1]);
            obj.setAcousticPaths('NS', NSacPathIR);

            obj.domain= 'frequency';
            NSacPathFR = repmat(permute(NS_FR(:, :, ns, rt), [1 3 2]), [1 2 1]);
            NSacPathFRstruct = struct('acousticPaths', NSacPathFR, 'frequencies', freqs);
            obj.setAcousticPaths('NS', NSacPathFRstruct); % obj.setAcousticPaths('NS', 'theoretical');

            if timeDomainActive
            % Simulate for time domain
            obj.domain = 'time';
            
                obj.WFSToolObj.filtersWFS_IR = repmat(permute(filtersWFS_IR(:, :, ns), [1 3 2]), [1 2 1]);
                            
                % Simulate only the noise source
                obj.WFSToolObj.virtual = [false; false];
                obj.WFSToolObj.WFScalculation();
                obj.WFSToolObj.simulate();
                recNS_signal = obj.WFSToolObj.simulField;
                
                % Simulate all together
                obj.WFSToolObj.virtual = [false; true];
                obj.WFSToolObj.WFScalculation();
                obj.WFSToolObj.simulate();
                rec_signal = obj.WFSToolObj.simulField;
                                
                % Identify IQ component
                recNScoef_time(:, f, ns, rt) = signal2pulseCoefficientMatrix([0 durSign], fcurr, 1, recNS_signal', fs).';
                recWFScoef_time(:, f, ns, rt) = signal2pulseCoefficientMatrix([0 durSign], fcurr, 1, (rec_signal - recNS_signal)', fs);
                recCoef_time(:, f, ns, rt) = signal2pulseCoefficientMatrix([0 durSign], fcurr, 1, rec_signal', fs);
                WFScoef_time(:, f, ns, rt) = signal2pulseCoefficientMatrix([0 durSign], fcurr, 1, obj.WFSToolObj.WFSarrayCoefficient', fs);
            end
                
            if frequencyDomainActive
            % Simulate for frequency domain
            obj.domain = 'frequency';

                obj.cancelResults = [];

                % Simulate 
                obj.WFSToolObj.virtual = [false; true];
                obj.cancel();

                recNScoef_freq(:, f, ns, rt) = obj.cancelResults.recNScoef;
                recCoef_freq(:, f, ns, rt) = obj.cancelResults.recCoef;
                recWFScoef_freq(:, f, ns, rt) = obj.cancelResults.recWFScoef;
                WFScoef_freq(:, f, ns, rt) = obj.cancelResults.WFScoef;
            end
            
                % Progress monitoring
                completedIterPerLevel = [rt-1, f-1, ns]; % Completed
                for k = numLevels:-1:2
                    if completedIterPerLevel(k) == maxIterPerLevel(k)
                        completedIterPerLevel(k) = 0;
                        completedIterPerLevel(k - 1) = completedIterPerLevel(k - 1) + 1;
                    end
                end
                completedIter = completedIterPerLevel*levelWeight;
                progress = completedIter/totalIter;
                
                aux = [completedIterPerLevel; maxIterPerLevel];
                msg = ['[', sprintf(' %d/%d ', aux), ']'];
                
                waitbar(progress, wb, msg);
                
                
        end
        
    end
end


%% Formatting of results

% We introduced an initial delay equal to the frequency filter delay, so we
% must correct the phase of the time coefficients.
phaseShift = filterDelay/fs * 2*pi*freqs;
recNScoef_time_correct = recNScoef_time .* exp(1i* repmat(phaseShift, [obj.numMicro, 1, numNSpos, numReverbTime]) );
recWFScoef_time_correct = recWFScoef_time .* exp(1i* repmat(phaseShift, [obj.numMicro, 1, numNSpos, numReverbTime]) );
recCoef_time_correct = recCoef_time .* exp(1i* repmat(phaseShift, [obj.numMicro, 1, numNSpos, numReverbTime]) );
WFScoef_time_correct = WFScoef_time .* exp(1i* repmat(phaseShift, [obj.numWFS, 1, numNSpos, numReverbTime]) );

% Make structure so we can visualize this
s_time = repmat(obj.cancelResults(1), [1, numFreqs, numNSpos, numReverbTime]);
s_freq = repmat(obj.cancelResults(1), [1, numFreqs, numNSpos, numReverbTime]);
for rt = 1:numReverbTime
    for f = 1:numFreqs
        for ns = 1:numNSpos
            s_time(1, f, ns, rt) = SimulationController.generateExportStructure(...
                'NSRcoef', 1,...
                'NSVcoef', -1,...
                'WFScoef', WFScoef_time_correct(:, f, ns, rt),...
                'microCoef', recCoef_time_correct(:, f, ns, rt),...
                'microCoefNS', recNScoef_time_correct(:, f, ns, rt),...
                'microCoefWFS', recWFScoef_time_correct(:, f, ns, rt),...
                'NSRpos', NSpositions(ns,:),...
                'NSVpos', NSpositions(ns,:),...
                'WFSpos', obj.WFSposition,...
                'microPos', obj.microPos,...
                'Frequency', freqs(f)...
                );
            
            s_freq(1, f, ns, rt) = SimulationController.generateExportStructure(...
                'NSRcoef', 1,...
                'NSVcoef', -1,...
                'WFScoef', WFScoef_freq(:, f, ns, rt),...
                'microCoef', recCoef_freq(:, f, ns, rt),...
                'microCoefNS', recNScoef_freq(:, f, ns, rt),...
                'microCoefWFS', recWFScoef_freq(:, f, ns, rt),...
                'NSRpos', NSpositions(ns,:),...
                'NSVpos', NSpositions(ns,:),...
                'WFSpos', obj.WFSposition,...
                'microPos', obj.microPos,...
                'Frequency', freqs(f)...
                );
        end
    end
end

s = s_freq; 

% Format structure
for p = 1:numel(s)
    s(p).NScoef = [s(p).NSRcoef; s(p).NSVcoef];
    s(p).NSposition = [s(p).NSRposition; s(p).NSVposition];
end

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

% Represent
% Graph 1
% ax = axes(figure);
% ax.Title.Interpreter = 'latex';
% hold on
% data = s(1).corrFactIndividual;
% cmap = colormap('lines');
% for k = 1:size(data, 2)
% xData = real(data(:, k));
% yData = imag(data(:, k));
% scat = scatter(ax, xData, yData, 10, cmap(k, :), 'filled');
% end
% maxAbs = max(abs(data(:)));
% ax.XLim = [-maxAbs, maxAbs]; ax.YLim = [-maxAbs, maxAbs];
% ax.DataAspectRatio = [1 1 3];
% 
% % Remember that mean and std are performed over the first dimension
% meanData = mean(data);
% stdData = std(data);
% normStd = stdData./abs(meanData);
% 
% meanAbs = mean(abs(data));
% stdAbs = std(abs(data));
% normStdAbs = stdAbs./meanAbs;
% 
% [meanPhase, ~, stdPhase] = circularDistributionParameters(angle(data));
% stdPhase = rad2deg(stdPhase);
% meanPhase = rad2deg(meanPhase);
% normStdPhase = stdPhase./meanPhase;

% Check the evolution of the individual correction factor with frequency.
% It should be similar for all points
[ ~, corrFactIndFirstFreq ] = pointWiseExtend( corrFactInd, corrFactInd(:, 1, :, :, :) );
relCorrFactInd = corrFactInd./corrFactIndFirstFreq;

% % Graph 2
% % Visualize. Don't delete just in case I need to use at another time.
% domain = [0];
% visualObj = animation({domain, freqs, 1:numNSpos, 1:numReverbTime, 1:numMicro},...
%     {abs(relCorrFactInd)}, {'Domain', 'Frequency', 'NS position index', 'Reverb. Time index', 'Micro Index'}, ...
%     {'Relative individual correction factor'}, [], []);
% 
% visualObj = animation({domain, freqs, 1:numNSpos, 1:numReverbTime, 1:numMicro},...
%     {abs(corrFactInd)}, {'Domain', 'Frequency', 'NS position index', 'Reverb. Time index', 'Micro Index'}, ...
%     {'Individual correction factor'}, [], []);
% 
% visualObj = animation({domain, freqs, 1:numNSpos, 1:numReverbTime},...
%     {abs(corrFactGlobal)}, {'Domain', 'Frequency', 'NS position index', 'Reverb. Time index'}, {'Cancellation'}, [], []);

% Graph 3
% How disperse are the individual correction factors for different NS
% positions and microphone positions?
% Note: you can use histogram2 instead of pcolor for visualizing, but the
% flexibility is smaller
freqs_aux = 0:1000;
correctFactTheo = sqrt(1i * freqs_aux/c);

corrFactInd_zeroReverb = corrFactInd(1, :, :, 1, :);
corrFactGlobal_zeroReverb = corrFactGlobal(1, :, :, 1);
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

Cg_theo_dB_FS = Cg_theo_dB(1, :, :, 1); % FS: Free Space
Cg_opt_dB_FS = Cg_opt_dB(1, :, :, 1); % FS: Free Space

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

Cg_A_dB_FS = Cg_A_dB(1, :, :, 1); % FS: Free Space
Cg_B_dB_FS = Cg_B_dB(1, :, :, 1); % FS: Free Space

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
magnFiltOrder = 512;
hilbertFiltOrder = 1024;
numFreqFilters = length(magnFiltOrder);

freqFilters = cell(numFreqFilters, 1);
freqFiltDelays = zeros(numFreqFilters, 1);
freqFiltFR = zeros(numFreqFilters, numFreqs);
for k = 1:numFreqFilters
    [freqFilter, delay] = getFrequencyFilter( magnFiltOrder(k), hilbertFiltOrder(k), fs );
    freqFiltFR(k, :) = DFT_slow(fs*freqFilter', fs, freqs).';
    
    freqFilters{k} = freqFilter;
    freqFiltDelays(k) = delay;
end
phaseShifts = -freqFiltDelays/fs * 2*pi*freqs;

% Simulation
durSign = 1; % Duration of signal
t = (0:ceil(durSign*obj.Fs)-1)/obj.Fs;
signalLength = length(t);
[~, filterDelay] = max(obj.WFSToolObj.freqFilter);
filterDelay = filterDelay - 1;

recNScoef_time = zeros(numMicro, numFreqs, numNSpos, numReverbTime, numFreqFilters);
recCoef_time = zeros(numMicro, numFreqs, numNSpos, numReverbTime, numFreqFilters);
recWFScoef_time = zeros(numMicro, numFreqs, numNSpos, numReverbTime, numFreqFilters);
WFScoef_time = zeros(obj.numWFS, numFreqs, numNSpos, numReverbTime, numFreqFilters);

recNScoef_freq = zeros(obj.numMicro, numFreqs, numNSpos, numReverbTime);
recCoef_freq = zeros(obj.numMicro, numFreqs, numNSpos, numReverbTime);
recWFScoef_freq = zeros(obj.numMicro, numFreqs, numNSpos, numReverbTime);
WFScoef_freq = zeros(obj.numWFS, numFreqs, numNSpos, numReverbTime);

timeDomainActive = true;
frequencyDomainActive = true;

obj.WFSToolObj.frequencyCorrection = true;

% Progress bar variables
maxIterPerLevel = [numReverbTime, numFreqs, numNSpos];
numLevels = length(maxIterPerLevel);
totalIter = prod(maxIterPerLevel);
levelWeight = flip([1, cumprod(flip(maxIterPerLevel(2:end)))])';
% completedIterPerLevel = [0 0 0]; % Completed
% completedIter = completedIterPerLevel*levelWeight;
% progress = completedIter/totalIter;
wb = waitbar(0, 'Progress');

for rt = 1:numReverbTime
    
    % Set up acoustic paths for the WFS array
    if WFS_AcPath_previously_calculated
        % In case you calculate WFS_IR and WFS_FR previously
        WFSacPathIR = permute(WFS_IR(:, :, :, rt), [1, 3, 2]);
        WFSacPathFR = permute(WFS_FR(:, :, :, rt), [1, 3, 2]);
    else
        WFS_IR_current = zeros(obj.numMicro, numSampIR, obj.numWFS);
        for k = 1:obj.numWFS
            WFS_IR_current(:, :, k) = rir_generator(c, fs, r, wfsPos(k, :), roomDim, Beta(rt, :), numSampIR);
        end
        WFSacPathIR = permute(WFS_IR_current, [1, 3, 2]);
        
        WFS_FR_current = zeros(obj.numMicro, numFreqs, obj.numWFS);
        for k = 1:obj.numWFS
            WFS_FR_current(:, :, k) = DFT_slow(fs*WFS_IR_current(:, :, k).', fs, freqs).';
        end
        WFSacPathFR = permute(WFS_FR_current, [1, 3, 2]);
    end
    
    obj.domain= 'time';
    obj.setAcousticPaths('WFS', WFSacPathIR);
    
    obj.domain= 'frequency';
    WFSacPathFRstruct = struct('acousticPaths', WFSacPathFR, 'frequencies', freqs);
    obj.setAcousticPaths('WFS', WFSacPathFRstruct); % obj.setAcousticPaths('WFS', 'theoretical');
    
    for f = 1:numFreqs
        fcurr = freqs(f); % Current frequency
        obj.frequency = fcurr;
        
        % Generate tone with the propper frequency
        x = obj.amplitude(1) * cos(2*pi*fcurr*t + obj.phase(1));
        x = [zeros(1, filterDelay), x, zeros(1, numSampIR - 1)];
               
        obj.domain = 'time';
        obj.NScoef = x;
        obj.NSVcoef = -x;
        
        for ns = 1:numNSpos         
            obj.NSposition = NSpositions(ns, :);
            
            % Set up acoustic paths for the noise source
            obj.domain= 'time';
            NSacPathIR = repmat(permute(NS_IR(:, :, ns, rt), [1, 3, 2]), [1, 2, 1]);
            obj.setAcousticPaths('NS', NSacPathIR);

            obj.domain= 'frequency';
            NSacPathFR = repmat(permute(NS_FR(:, :, ns, rt), [1 3 2]), [1 2 1]);
            NSacPathFRstruct = struct('acousticPaths', NSacPathFR, 'frequencies', freqs);
            obj.setAcousticPaths('NS', NSacPathFRstruct); % obj.setAcousticPaths('NS', 'theoretical');

            % Time domain calculation. FAKE TIME PROCESSING.
            if timeDomainActive
                % Simulate for time domain
%                 obj.domain = 'time';
%             
%                 obj.WFSToolObj.filtersWFS_IR = repmat(permute(filtersWFS_IR(:, :, ns), [1 3 2]), [1 2 1]);
%                             
%                 % Simulate only the noise source
%                 obj.WFSToolObj.virtual = [false; false];
%                 obj.WFSToolObj.WFScalculation();
%                 obj.WFSToolObj.simulate();
%                 recNS_signal = obj.WFSToolObj.simulField;
%                 
%                 % Simulate all together
%                 obj.WFSToolObj.virtual = [false; true];
%                 obj.WFSToolObj.WFScalculation();
%                 obj.WFSToolObj.simulate();
%                 rec_signal = obj.WFSToolObj.simulField;
%                 
%                 % Identify IQ component
%                 recNScoef_time(:, f, ns, rt) = signal2pulseCoefficientMatrix([0 durSign], fcurr, 1, recNS_signal', fs).';
%                 recWFScoef_time(:, f, ns, rt) = signal2pulseCoefficientMatrix([0 durSign], fcurr, 1, (rec_signal - recNS_signal)', fs);
%                 recCoef_time(:, f, ns, rt) = signal2pulseCoefficientMatrix([0 durSign], fcurr, 1, rec_signal', fs);
%                 WFScoef_time(:, f, ns, rt) = signal2pulseCoefficientMatrix([0 durSign], fcurr, 1, obj.WFSToolObj.WFSarrayCoefficient', fs);
                
                % We introduce a modification because previous code is
                % innefficient.
                % We know that the only difference between the processing
                % in frequency and time is, besides the fact that the time signals
                % are finite and hence there will be truncatin errors, that
                % the frequency response of the frequency filters for the
                % calculation of WFS signals is not exactly equal to the
                % ideal case.
                                
                obj.domain = 'frequency';
                obj.WFSToolObj.virtual = [false; true];
                
                obj.WFSToolObj.frequencyCorrection = false;
                obj.WFSToolObj.WFScalculation();
                obj.WFSToolObj.frequencyCorrection = true;
                
                for filt = 1:numFreqFilters
                    obj.WFSToolObj.WFSarrayCoefficient = obj.WFSToolObj.WFSarrayCoefficient * freqFilterFR(f) * exp(-1i*phaseShift(f));
                    obj.WFSToolObj.simulate();
                    
                    WFScoef_time(:, f, ns, rt, filt) = obj.WFScoef;
                    recCoef_time(:, f, ns, rt, filt) = obj.microCoef;
                    recWFScoef_time(:, f, ns, rt, filt) = obj.microCoefWFS;
                    recNScoef_time(:, f, ns, rt, filt) = obj.microCoefNS;
                end
                
                                
                
                
            end
                
            % Frequency domain calculation
            if frequencyDomainActive
            % Simulate for frequency domain
            obj.domain = 'frequency';

                obj.cancelResults = [];

                % Simulate 
                obj.WFSToolObj.virtual = [false; true];
                obj.cancel();

                recNScoef_freq(:, f, ns, rt) = obj.cancelResults.recNScoef;
                recCoef_freq(:, f, ns, rt) = obj.cancelResults.recCoef;
                recWFScoef_freq(:, f, ns, rt) = obj.cancelResults.recWFScoef;
                WFScoef_freq(:, f, ns, rt) = obj.cancelResults.WFScoef;
            end
            
                % Progress monitoring
                completedIterPerLevel = [rt-1, f-1, ns]; % Completed
                for k = numLevels:-1:2
                    if completedIterPerLevel(k) == maxIterPerLevel(k)
                        completedIterPerLevel(k) = 0;
                        completedIterPerLevel(k - 1) = completedIterPerLevel(k - 1) + 1;
                    end
                end
                completedIter = completedIterPerLevel*levelWeight;
                progress = completedIter/totalIter;
                
                aux = [completedIterPerLevel; maxIterPerLevel];
                msg = ['[', sprintf(' %d/%d ', aux), ']'];
                
                waitbar(progress, wb, msg);
                
                
        end
        
    end
end

%% Formatting of results

% Make structure so we can visualize this
s_time = repmat(obj.cancelResults(1), [1, numFreqs, numNSpos, numReverbTime]);
s_freq = repmat(obj.cancelResults(1), [1, numFreqs, numNSpos, numReverbTime]);
for rt = 1:numReverbTime
    for f = 1:numFreqs
        for ns = 1:numNSpos
            s_time(1, f, ns, rt) = SimulationController.generateExportStructure(...
                'NSRcoef', 1,...
                'NSVcoef', -1,...
                'WFScoef', WFScoef_time(:, f, ns, rt),...
                'microCoef', recCoef_time(:, f, ns, rt),...
                'microCoefNS', recNScoef_time(:, f, ns, rt),...
                'microCoefWFS', recWFScoef_time(:, f, ns, rt),...
                'NSRpos', NSpositions(ns,:),...
                'NSVpos', NSpositions(ns,:),...
                'WFSpos', obj.WFSposition,...
                'microPos', obj.microPos,...
                'Frequency', freqs(f)...
                );
            
            s_freq(1, f, ns, rt) = SimulationController.generateExportStructure(...
                'NSRcoef', 1,...
                'NSVcoef', -1,...
                'WFScoef', WFScoef_freq(:, f, ns, rt),...
                'microCoef', recCoef_freq(:, f, ns, rt),...
                'microCoefNS', recNScoef_freq(:, f, ns, rt),...
                'microCoefWFS', recWFScoef_freq(:, f, ns, rt),...
                'NSRpos', NSpositions(ns,:),...
                'NSVpos', NSpositions(ns,:),...
                'WFSpos', obj.WFSposition,...
                'microPos', obj.microPos,...
                'Frequency', freqs(f)...
                );
        end
    end
end

s = [s_time; s_freq]; 

% Format structure
for p = 1:numel(s)
    s(p).NScoef = [s(p).NSRcoef; s(p).NSVcoef];
    s(p).NSposition = [s(p).NSRposition; s(p).NSVposition];
end

%% Global cancellation with the filter response

Cg = zeros(size(s));
for k = 1:numel(s)
    Cg(k) = sum(abs(s(k).recCoef).^2)/sum(abs(s(k).recNScoef).^2);
end
Cg_dB = 10*log10(Cg);

% Visualize
freqMat_FS = pointWiseExtend(freqs, s(1, :, :, 1));

Cg_dB_FS_time = Cg_dB(1, :, :, 1); % Time domain. FS: Free Space
Cg_dB_FS_freq = Cg_dB(2, :, :, 1); % Frequency domain. FS: Free Space

freqStep = freqs(2) - freqs(1);
freqEdges = [(freqs(1:end-1) + freqs(2:end))/2, freqs(end) + freqStep/2];
CancEdges = -20:10;

N_time = histcounts2(freqMat_FS, Cg_dB_FS_time, freqEdges, CancEdges);
Nnorm_time = N_time./repmat(sum(N_time, 2), [1, size(N_time, 2)]); Nnorm_time(isnan(Nnorm_time)) = 0;

N_freq = histcounts2(freqMat_FS, Cg_dB_FS_freq, freqEdges, CancEdges);
Nnorm_freq = N_freq./repmat(sum(N_freq, 2), [1, size(N_freq, 2)]); Nnorm_freq(isnan(Nnorm_freq)) = 0;

penaltyEdges = -10:10;
[N_diff, freqAux, DifCancEdges] = histcounts2(freqMat_FS, Cg_dB_FS_time - Cg_dB_FS_freq, freqEdges, penaltyEdges);
Nnorm_diff = N_diff./repmat(sum(N_diff, 2), [1, size(N_diff, 2)]);


ax1 = axes(figure);
C = zeros(length(freqEdges), length(CancEdges));
C(1:end-1, 1:end-1) = Nnorm_time;
pcolor(ax1, freqEdges, CancEdges, C');
ax1.Title.String = ['$C_{global}$ for time processing'];
ax1.Title.Interpreter = 'latex';
ax1.XLabel.String = 'Frequency (Hz)';
ax1.YLabel.String = '$C_{global}$';
ax1.YLabel.Interpreter = 'latex';
colorbar(ax1)
% printfig(ax1.Parent, imagesPath, 'Experiment7_CancGlobTheoCorrFactTime', 'eps');

ax2 = axes(figure);
C = zeros(length(freqEdges), length(CancEdges));
C(1:end-1, 1:end-1) = Nnorm_freq;
pcolor(ax2, freqEdges, CancEdges, C');
ax2.Title.String = ['$C_{global}$ for frequency processing'];
ax2.Title.Interpreter = 'latex';
ax2.XLabel.String = 'Frequency (Hz)';
ax2.YLabel.String = '$C_{global}$';
ax2.YLabel.Interpreter = 'latex';
colorbar(ax2)
% printfig(ax2.Parent, imagesPath, 'Experiment7_CancGlobTheoCorrFactFreq', 'eps');

ax3 = axes(figure);
C = zeros(length(freqEdges), length(penaltyEdges));
C(1:end-1, 1:end-1) = Nnorm_diff;
pcolor(ax3, freqEdges, penaltyEdges, C');
ax3.Title.String = ['Cancellation Penalty'];
ax3.Title.Interpreter = 'latex';
ax3.XLabel.String = 'Frequency (Hz)';
ax3.YLabel.String = '$10\log_{10}(C_{global,time}) - 10\log_{10}(C_{global,ideal})$';
ax3.YLabel.Interpreter = 'latex';
colorbar(ax3)


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


