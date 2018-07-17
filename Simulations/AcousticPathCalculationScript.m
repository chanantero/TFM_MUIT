%% Pre-calculate impulse responses
% Since it would take a lot of time to calculate the IR during the
% simulations, we calculate it previously and save it to a .mat. 

r = recPositions + repmat(WFSarrayOffset, numMicro, 1); % Receiver position [x y z] (m)
wfsPos = obj.WFSposition + repmat(WFSarrayOffset, obj.numWFS, 1);
nsPos = NSpositions + repmat(WFSarrayOffset, numNSpos, 1);

if ~predefRoomDim
    maxX = max([nsPos(:, 1); wfsPos(:, 1)]);
    maxY = max([nsPos(:, 2); wfsPos(:, 2)]);
    roomDim = [maxX+1 maxY+1 4];                % Room dimensions [x y z] (m)
end

if ~predefNumSampIR
    % Minimize the number of samples
    dist = zeros(numMicro, numNSpos);
    for ns = 1:numNSpos
        dist(:, ns) = sqrt(sum((recPositions - repmat(NSpositions(ns, :), [numMicro, 1])).^2, 2));
    end
    maxDist = max(dist(:));
    numSampIR = 2^ceil(log2(maxDist/c*fs)); % Number of samples
end

if ~predefWFSfilterLength
    WFSfilterLength = numSampIR;
    obj.WFSToolObj.filterWFS_length = WFSfilterLength;
end

if WFS_AcPath_previously_calculated
    WFS_IR = zeros(obj.numMicro, numSampIR, obj.numWFS, numReverbTime);
    for k = 1:obj.numWFS
        for rt = 1:numReverbTime
            WFS_IR(:, :, k, rt) = rir_generator(c, fs, r, wfsPos(k, :), roomDim, Beta(rt, :), numSampIR);
        end
    end
    
    % Calculate the frequency responses
%     WFS_FR = zeros(obj.numMicro, numFreqs, obj.numWFS, numReverbTime);
%     for k = 1:obj.numWFS
%         for rt = 1:numReverbTime
%             WFS_FR(:, :, k, rt) = DFT_slow(fs*WFS_IR(:, :, k, rt).', fs, freqs).';
%         end
%     end
    
    if numFreqs > 1
        oper = @(x) freqz(x, 1, freqs, fs);
        WFS_FR = oneDimOperOverMultiDimArray( oper, WFS_IR, 2 );
    else
        oper = @(x) freqz(x, 1, [0 freqs], fs);
        WFS_FR = oneDimOperOverMultiDimArray( oper, WFS_IR, 2 );
        WFS_FR = WFS_FR(:, 2, :, :);
    end
    
end

% Do the same with the noise source positions
if NS_AcPath_previously_calculated
    NS_IR = zeros(obj.numMicro, numSampIR, numNSpos, numReverbTime);
    for k = 1:numNSpos
        for rt = 1:numReverbTime
            NS_IR(:, :, k, rt) = rir_generator(c, fs, r, nsPos(k,:), roomDim, Beta(rt, :), numSampIR);
        end
    end
    
%     NS_FR = zeros(obj.numMicro, numFreqs, numNSpos, numReverbTime);
%     for k = 1:numNSpos
%         for rt = 1:numReverbTime
%             NS_FR(:, :, k, rt) = DFT_slow(fs*NS_IR(:, :, k, rt).', fs, freqs).';
%         end
%     end
    
    oper = @(x) freqz(x, 1, freqs, fs);
    NS_FR = oneDimOperOverMultiDimArray( oper, NS_IR, 2 );
    
end

if WFS_AcPath_previously_calculated && NS_AcPath_previously_calculated ...
        && appendFreeSpaceAcPaths
    
    if timeDomainActive
        acPath = simulator.calculateMonopolesIR(obj.WFSposition, recPositions, c, fs, numSampIR);
        WFS_IR = cat(4, permute(acPath, [1 3 2 4]), WFS_IR);
    end
    
    acPath = simulator.calculateTheoricAcousticPaths(...
        obj.WFSToolObj.WFSarrayPosition, obj.WFSToolObj.WFSarrayRadiationPattern, obj.WFSToolObj.WFSarrayOrientation,...
        recPositions, obj.WFSToolObj.receiverRadiationPattern, obj.WFSToolObj.receiverOrientation, freqs, c);
    WFS_FR = cat(4, permute(acPath, [1 3 2 4]), WFS_FR);
    
    if timeDomainActive
        acPath = simulator.calculateMonopolesIR(NSpositions, recPositions, c, fs, numSampIR);
        NS_IR = cat(4, permute(acPath, [1 3 2 4]), NS_IR);
    end
    
    acPath = simulator.calculateTheoricAcousticPaths(...
        NSpositions, repmat({@simulator.monopoleRadPat}, numNSpos, 1), repmat([1 0 0 1], numNSpos, 1),... % [1 0 0 1] is just a default orientation. It doesn't matter because they are monopoles
        recPositions, obj.WFSToolObj.receiverRadiationPattern, obj.WFSToolObj.receiverOrientation, freqs, c);
    NS_FR = cat(4, permute(acPath, [1 3 2 4]), NS_FR);
    
    numReverbTime = numReverbTime + 1;
end