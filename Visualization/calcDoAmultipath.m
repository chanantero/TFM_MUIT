function [DoAs_ULA, DoAs_sp, DoAs_ne, DoAs_cs, varargout] = calcDoAmultipath( transDoAs, signalPowers, multipath, elecDist, numAntennas, processingEnvironment )
% Everything ideal except multipath.
% Inputs:
% - transDoAs. Horizontal vector with the Directions of Arrival of
% transmitters.
% - signalAmplitudes. Horizontal vector with the amplitudes of the signals.
% - multipath. Horizontal cell vector with as many elements as sources.
% The i-th cell contains a matrix with 3 rows and as many columns as multipath
% components has the i-th source. The first row contains the DoA of the
% multipath components, the second one the relative amplitude with respect
% to the real signal, and the third one the relative phase. The phases and
% amplitudes are relative to the center of the receiver ULA.
% - elecDist. Separation of the antennas of the receiver Uniform Linear Array
% in wavelengths.
% - numAntennas. Number of receiver antennas.
% - processingEnvironment.

outOfBounds = false;

numSources = numel(signalPowers); % numSources = numel(transDoAs);

numMultipath = cellfun(@(mult) size(mult,2), multipath);
numVirtualSources = sum(numMultipath);
numTotalSources = numSources + numVirtualSources;
A = zeros(numAntennas, numTotalSources);
Rss = zeros(numTotalSources);
last = 0;
for s = 1:numSources
    multipathMat = multipath{s};

    relAmplitudes = [1; multipathMat(2,:)'];
    relPhases = [0; multipathMat(3,:).'];
    % delay of the multipath signals respect the LOS one, being 1 the symbol period
    if(size(multipathMat, 1) < 4)
        delays = zeros(1+numMultipath(s), 1);
    else
        delays = [0; multipathMat(4,:)'];
    end
    
    relAmpMat = relAmplitudes*relAmplitudes';
    relPhaMat = exp(1i*relPhases)*exp(1i*relPhases)';
    [a, b] = ndgrid(delays, delays);
    corrRatMat = max(1 - abs(a - b), 0);
    
    Rss(last+1:last+1+numMultipath(s), last+1:last+1+numMultipath(s)) = signalPowers(s)*corrRatMat.*relAmpMat.*relPhaMat;
    
    mainStVec = steeringVector(numAntennas, elecDist, transDoAs(s), true);
    multStVec = zeros(numAntennas, numMultipath(s));

    for vs = 1:numMultipath(s)
        multDoA = multipathMat(1,vs);
        if(any(abs(multDoA) > 90)) % Out of [-90°, 90°] bounds
            outOfBounds = true;
            break;
        end
        multStVec(:,vs) = steeringVector(numAntennas, elecDist, multDoA, true);
    end
    if(outOfBounds);  break;  end
    
    A(:, last+1:last+(1+numMultipath(s))) = [mainStVec, multStVec];
    last = last + (1+numMultipath(s));
end

if(outOfBounds)
    DoAs_ULA = [];
    DoAs_sp = [];
    DoAs_ne = [];
    DoAs_cs = [];
    return;
end

% R3 = A*Rss*A';
% En vez de hacer una multiplicación de matrices habitual, usar un bucle
% for
R = zeros(numAntennas, numAntennas, numTotalSources^2);
for row = 1:numTotalSources
    for column = 1:numTotalSources
        R(:, :, (row-1)*numTotalSources + column) = A(:,row)*A(:,column)'*Rss(row, column);
    end
end
R = sum(R, 3);

% % Otras formas que no llegaron a funcionar bien, o bien por falta de
% funcionalidad o bien porque la matriz resultante era casi singular y daba
% problema al hacer la proyección en el eigenspace de ruido
% numMultipath = cellfun(@(mult) size(mult,2), multipath);
% numVirtualSources = sum(numMultipath);
% numTotalSources = numSources + numVirtualSources;
% A = zeros(numAntennas, numTotalSources);
% Rss = zeros(numTotalSources);
% last = 0;
% for s = 1:numSources
%     multipathMat = multipath{s};
% 
%     relAmplitudes = multipathMat(2,:)';
%     relPhases = multipathMat(3,:).';
%     % delay of the multipath signals respect the LOS one, being 1 the symbol period
%     if(size(multipathMat, 1) < 4)
%         delays = zeros(1+numMultipath(s), 1);
%     else
%         delays = [0; multipathMat(4,:)'];
%     end
%     
%     [a, b] = ndgrid(delays, delays);
%     corrRatMat = 1 - abs(a - b);
%     
%     Rss(last+1:last+1+numMultipath(s), last+1:last+1+numMultipath(s)) = signalPowers(s)*corrRatMat;
%     
%     mainStVec = steeringVector(numAntennas, elecDist, transDoAs(s), true);
%     multStVec = zeros(numAntennas, numMultipath(s));
% 
%     for vs = 1:numMultipath(s)
%         multDoA = multipathMat(1,vs);
%         if(any(abs(multDoA) > 90)) % Out of [-90°, 90°] bounds
%             outOfBounds = true;
%             break;
%         end
%         multStVec(:,vs) = steeringVector(numAntennas, elecDist, multDoA, true)*relAmplitudes(vs)*exp(1i*relPhases(vs));
%     end
%     if(outOfBounds);  break;  end
%     
%     A(:, last+1:last+(1+numMultipath(s))) = [mainStVec, multStVec];
%     last = last + (1+numMultipath(s));
% end
% 
% if(outOfBounds)
%     DoAs_ULA = [];
%     DoAs_sp = [];
%     DoAs_ne = [];
%     DoAs_cs = [];
%     return;
% end
% 
% R2 = A*Rss*A';
% 
% % Forma simple: no tiene en cuenta la posibilidad de que la correlación
% % entre señales redtardadas no sea compleata
% 
% % For each source, we calculate the steering vector of the LOS signal and
% % the multipath components. Then we sum them accordingly and include them
% % in the array manifold
% A = zeros(numAntennas, numSources); % Array manifold
% for s = 1:numSources
%     mainStVec = steeringVector(numAntennas, elecDist, transDoAs(s), true);
%     
%     multipathMat = multipath{s};
%     numMultipath = size(multipathMat, 2);
%     multStVec = zeros(numAntennas, numMultipath);
%     for m = 1:numMultipath
%         multDoA = multipathMat(1,m);
%         if(any(abs(multDoA) > 90)) % Out of [-90°, 90°] bounds
%             outOfBounds = true;
%             break;
%         end
%         relAmplitude = multipathMat(2,m);
%         relPhase = multipathMat(3,m);
%         multStVec(:,m) = steeringVector(numAntennas, elecDist, multDoA, true)*relAmplitude*exp(1i*relPhase);
%     end
%     if(outOfBounds);  break;  end
%     
%     A(:,s) = mainStVec + sum(multStVec, 2);
% end
% 
% if(outOfBounds)
%     DoAs_ULA = [];
%     DoAs_sp = [];
%     DoAs_ne = [];
%     DoAs_cs = [];
%     return;
% end
% 
% % Correlation matrix
% R = zeros(numAntennas, numAntennas, numSources);
% for k = 1:numSources
%     R(:,:,k) = A(:,k)*A(:,k)'*signalPowers(k);
% end
% R1 = sum(R, 3);

measureData.Rxx = R;
measureData.L = numSources;
measureData.Delta = elecDist;

spectrumData = estimateDoA( measureData, processingEnvironment );

DoAs_ULA = [];
DoAs_sp = [];
DoAs_ne = [];
DoAs_cs = [];
DoAs_ULA_Height = [];
DoAs_sp_Height = [];
DoAs_ne_Height = [];
DoAs_cs_Height = [];
members = ismember({'disable_ULA', 'sparse_ruler', 'nested_elements', 'coprime_sampling'}, fields(processingEnvironment));
if(~members(1) || ~processingEnvironment.disable_ULA)
    DoAs_ULA = spectrumData.DoAs;
    DoAs_ULA_Height = spectrumData.DoAsHeight;
end
if(members(2))
    DoAs_sp = spectrumData.DoAs_sp;
    DoAs_sp_Height = spectrumData.DoAsHeight_sp;
end
if(members(3))
    DoAs_ne = spectrumData.DoAs_ne;
    DoAs_ne_Height = spectrumData.DoAsHeight_ne;
end
if(members(4))
    DoAs_cs = spectrumData.DoAs_cs;
    DoAs_cs_Height = spectrumData.DoAsHeight_cs;
end
if(nargout == 8)
    varargout = {DoAs_ULA_Height, DoAs_sp_Height, DoAs_ne_Height, DoAs_cs_Height};
end
if nargout == 9
    varargout = {DoAs_ULA_Height, DoAs_sp_Height, DoAs_ne_Height, DoAs_cs_Height, spectrumData};
end
end

