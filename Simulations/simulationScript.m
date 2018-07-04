%% Simulation
% Creation of frequency filters with different orders
freqFiltFR = zeros(numFreqFilters, numFreqs);
for k = 1:numFreqFilters
%     freqFiltFR(k, :) = DFT_slow(fs*freqFilters{k}', fs, freqs).';
    freqFiltFR(k, :) = freqz(freqFilters{k}, 1, freqs, fs);
end
phaseShifts = -freqFiltDelays/fs * 2*pi*freqs;

% Signals
if timeDomainActive && ~fakeTimeProcessing
    if ~predefSignals
        t = (0:ceil(durSign*obj.Fs)-1)/obj.Fs;
        x = obj.amplitude(1) * cos(2*pi*freqs(:)*t + obj.phase(1));
    else
        x = NSsignal;
    end
    
    preDelay = max(freqFiltDelays);
    postDelay = numSampIR - 1 + max(freqFiltDelays);
    x = [zeros(numFreqs, preDelay), x, zeros(numFreqs, postDelay)];
    numSampSign = size(x, 2);
    
    if saveSignals
        recNS_signals = zeros(numMicro, numSampSign, numFreqs, numNSpos, numReverbTime);
        rec_signals = zeros(numMicro, numSampSign, numFreqs, numNSpos, numReverbTime, numFreqFilters);
    end
end

recNScoef_time = zeros(numMicro, numFreqs, numNSpos, numReverbTime);
recCoef_time = zeros(numMicro, numFreqs, numNSpos, numReverbTime, numFreqFilters);
recWFScoef_time = zeros(numMicro, numFreqs, numNSpos, numReverbTime, numFreqFilters);
WFScoef_time = zeros(obj.numWFS, numFreqs, numNSpos, numReverbTime, numFreqFilters);

recNScoef_freq = zeros(numMicro, numFreqs, numNSpos, numReverbTime);
recCoef_freq = zeros(numMicro, numFreqs, numNSpos, numReverbTime);
recWFScoef_freq = zeros(numMicro, numFreqs, numNSpos, numReverbTime);
WFScoef_freq = zeros(obj.numWFS, numFreqs, numNSpos, numReverbTime);

% Progress bar variables
if progressBarActive
maxIterPerLevel = [numReverbTime, numFreqs, numNSpos];
numLevels = length(maxIterPerLevel);
totalIter = prod(maxIterPerLevel);
levelWeight = flip([1, cumprod(flip(maxIterPerLevel(2:end)))])';
wb = waitbar(0, 'Progress');
descMsg = '';
aux = [zeros(1, numLevels); maxIterPerLevel];
numMsg = ['[', sprintf(' %d/%d ', aux), ']'];
progress = 0;
end
for rt = 1:numReverbTime
    
    % Set up acoustic paths for the WFS array
    if WFS_AcPath_previously_calculated
        % In case you calculate WFS_IR and WFS_FR previously
        if timeDomainActive
            WFSacPathIR = permute(WFS_IR(:, :, :, rt), [1, 3, 2]);
        end
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
    
    if timeDomainActive
    obj.domain= 'time';
    obj.setAcousticPaths('WFS', WFSacPathIR);
    end
    
    obj.domain= 'frequency';
    WFSacPathFRstruct = struct('acousticPaths', WFSacPathFR, 'frequencies', freqs);
    obj.setAcousticPaths('WFS', WFSacPathFRstruct); % obj.setAcousticPaths('WFS', 'theoretical');
    
    for f = 1:numFreqs
        fcurr = freqs(f); % Current frequency
        obj.frequency = fcurr;
        
        if timeDomainActive && ~fakeTimeProcessing
            preDelayPhaseShift = preDelay/fs*2*pi*fcurr;
        end
        
        % Generate tone with the propper frequency
%         x = obj.amplitude(1) * cos(2*pi*fcurr*t + obj.phase(1));
%         x = [zeros(1, preDelay), x, zeros(1, postDelay)];
        
        if timeDomainActive && ~fakeTimeProcessing
            obj.domain = 'time';
            obj.NScoef = x(f, :);
            obj.NSVcoef = -x(f, :);
        end
        
        for ns = 1:numNSpos
            obj.NSposition = NSpositions(ns, :);
            
            % Set up acoustic paths for the noise source
            if timeDomainActive
                obj.domain= 'time';
                NSacPathIR = repmat(permute(NS_IR(:, :, ns, rt), [1, 3, 2]), [1, 2, 1]);
                obj.setAcousticPaths('NS', NSacPathIR);
            end
            
            obj.domain= 'frequency';
            NSacPathFR = repmat(permute(NS_FR(:, :, ns, rt), [1 3 2]), [1 2 1]);
            NSacPathFRstruct = struct('acousticPaths', NSacPathFR, 'frequencies', freqs);
            obj.setAcousticPaths('NS', NSacPathFRstruct); % obj.setAcousticPaths('NS', 'theoretical');
            
            if timeDomainActive
%                 
                if ~fakeTimeProcessing
                    % Simulate for time domain
                    obj.domain = 'time';
                                        
                    % Simulate only the noise source
                    obj.WFSToolObj.virtual = [false; false];
                    obj.WFSToolObj.freqFilter = 1; % This line (¿and next line?) is necessary for technical issues. ¿Is it?
%                     obj.WFSToolObj.updateFiltersWFS();
                    obj.WFSToolObj.WFScalculation();
                    obj.WFSToolObj.simulate();
                    recNS_signal = obj.WFSToolObj.simulField;
                    recNScoef_time(:, f, ns, rt) = exp(1i*preDelayPhaseShift) * signal2pulseCoefficientMatrix([0 durSign], fcurr, 1, recNS_signal', fs).';
                    
                    if saveSignals
                        recNS_signals(:, :, f, ns, rt) = recNS_signal;
                    end                 
                    
                    obj.WFSToolObj.virtual = [false; true];
                    for filt = 1:numFreqFilters
                        % Simulate all together
                        obj.WFSToolObj.freqFilter = freqFilters{filt};
                        obj.WFSToolObj.updateFiltersWFS();
                        obj.WFSToolObj.WFScalculation();
                        obj.WFSToolObj.simulate();
                        rec_signal = obj.WFSToolObj.simulField;
                        
                        % Identify IQ component
                        recWFScoef_time(:, f, ns, rt, filt) = exp(1i*preDelayPhaseShift) * signal2pulseCoefficientMatrix([0 durSign], fcurr, 1, (rec_signal - recNS_signal)', fs);
                        recCoef_time(:, f, ns, rt, filt) = exp(1i*preDelayPhaseShift) * signal2pulseCoefficientMatrix([0 durSign], fcurr, 1, rec_signal', fs);
                        WFScoef_time(:, f, ns, rt, filt) = exp(1i*preDelayPhaseShift) * signal2pulseCoefficientMatrix([0 durSign], fcurr, 1, obj.WFSToolObj.WFSarrayCoefficient', fs);
                        
                        if saveSignals
                            rec_signals(:, :, f, ns, rt, filt) = rec_signal;
                        end
                        
                    end
                else
                    % Simulate for time domain
                    
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
                    WFScoef = obj.WFSToolObj.WFSarrayCoefficient;
                    obj.WFSToolObj.frequencyCorrection = true;
                    
                    for filt = 1:numFreqFilters
                        obj.WFSToolObj.WFSarrayCoefficient = WFScoef * freqFiltFR(filt, f) * exp(-1i*phaseShifts(filt, f));
                        obj.WFSToolObj.simulate();
                        
                        WFScoef_time(:, f, ns, rt, filt) = obj.WFScoef;
                        recCoef_time(:, f, ns, rt, filt) = obj.microCoef;
                        recWFScoef_time(:, f, ns, rt, filt) = obj.microCoefWFS;
                        recNScoef_time(:, f, ns, rt, filt) = obj.microCoefNS;
                    end
                end
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
            if progressBarActive % ns == numNSpos
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
                numMsg = ['[', sprintf(' %d/%d ', aux), ']'];
                msg = [numMsg, descMsg];
                
                waitbar(progress, wb, msg);
            end
            
        end
        
    end
end

if progressBarActive
close(wb);
end

%% Formatting of results

obj.domain = 'frequency';
sbase = obj.generateBasicExportStructure();

% Make structure so we can visualize this
s_time = repmat(sbase, [numFreqFilters, numFreqs, numNSpos, numReverbTime]);
s_freq = repmat(sbase, [1, numFreqs, numNSpos, numReverbTime]);
for rt = 1:numReverbTime
    for f = 1:numFreqs
        for ns = 1:numNSpos
            
            if timeDomainActive
                for filt = 1:numFreqFilters
                    s_time(filt, f, ns, rt) = SimulationController.generateExportStructure(...
                        'NSRcoef', 1,...
                        'NSVcoef', -1,...
                        'WFScoef', WFScoef_time(:, f, ns, rt, filt),...
                        'microCoef', recCoef_time(:, f, ns, rt, filt),...
                        'microCoefNS', recNScoef_time(:, f, ns, rt),...
                        'microCoefWFS', recWFScoef_time(:, f, ns, rt, filt),...
                        'NSRpos', NSpositions(ns,:),...
                        'NSVpos', NSpositions(ns,:),...
                        'WFSpos', obj.WFSposition,...
                        'microPos', obj.microPos,...
                        'Frequency', freqs(f)...
                        );
                end
            end
            
            if frequencyDomainActive
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
end 

if frequencyDomainActive && timeDomainActive
    s = [s_time; s_freq];
elseif frequencyDomainActive
    s = s_freq;
elseif timeDomainActive
    s = s_time;
end

% Format structure
for p = 1:numel(s)
    s(p).NScoef = [s(p).NSRcoef; s(p).NSVcoef];
    s(p).NSposition = [s(p).NSRposition; s(p).NSVposition];
end
