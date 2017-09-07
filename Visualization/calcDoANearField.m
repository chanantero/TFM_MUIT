function [DoAs_ULA, DoAs_sp, DoAs_ne, DoAs_cs, varargout] = calcDoANearField( transPosit, signalPowers, receivPosit, processingEnvironment )
% Everything ideal except the farfield. Propagation losses neglected
% The positions and distances are electric distances

% Transmitters
numSources = numel(signalPowers);
[~, ~, distancias] = sphericalCoord(transPosit, receivPosit);

% Receivers
numAntennas = size(receivPosit, 1);

% Correlation matrix
R = zeros(numAntennas, numAntennas, numSources);
for k = 1:numSources
    dist = distancias(:,k);
    diffDist = -repmat(dist, 1, numAntennas) + repmat(dist', numAntennas, 1);
    R(:,:,k) = exp(1i*2*pi*diffDist)*signalPowers(k);
end
R = sum(R, 3);

measureData.Rxx = R;
measureData.L = numSources;
measureData.Delta = 0.5;

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