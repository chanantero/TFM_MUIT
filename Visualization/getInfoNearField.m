function varargout = getInfoNearField(inputMethod, parameters, selectedOutput, varargin) % Varargin are extra output options
% getInfoMultipath( selectedSource, selectedOutput, selectedAlgorithm, maxDev, option, varargin )
% varargout = getInfoMultipath( mainDoA, signalPowers, multipath, elecDist, processingEnvironment, maxDev, selectedSource, selectedOutput, selectedAlgorithm )
% getInfoMultipath(inputMethod, parameters, outputOptions)

% Give form to inputs
[DoA, transPosit, signalPowers, receivPosit, processingEnvironment, selectedAlgorithm] = ...
    formatInputsNearField(inputMethod, parameters, 'extended');

% Calculate
[calcDoA_ULA, calcDoA_sp, calcDoA_ne, calcDoA_cs,...
 calcHeight_ULA, calcHeight_sp, calcHeight_ne, calcHeight_cs, spectrumData]...
 = calcDoANearField( transPosit, signalPowers, receivPosit, processingEnvironment );


investigated_theta = spectrumData.investigated_theta;
% Choose which algorithm to use
switch selectedAlgorithm
    case 'ULA'
        calcDoA = calcDoA_ULA;
        calcHeight = calcHeight_ULA;
        spectrum = spectrumData.norm_sp_mu;
    case 'SRA'
        calcDoA = calcDoA_sp;
        calcHeight = calcHeight_sp;
        spectrum = spectrumData.output_norm_sr;
    case 'NE'
        calcDoA = calcDoA_ne;
        calcHeight = calcHeight_ne;
        spectrum = spectrumData.output_norm_ne;
    case 'CS'
        calcDoA = calcDoA_cs;
        calcHeight = calcHeight_cs;
        spectrum = spectrumData.output_norm_cs;
end

if(strcmp(selectedOutput, 'spectrum'))
    spectrumData = struct;
    spectrumData.investigated_theta = investigated_theta;
    spectrumData.spectrum = spectrum;
    spectrumData.calcDoA = calcDoA;
    varargout{1} = spectrumData;
else
    varargout{1} = getPeaksInfo(DoA, calcDoA, calcHeight, {selectedOutput}, varargin);
end

end

% function varargout = formatInputs(option, parameters)
% switch option
%     case '1source'
%         [DoA, dist, selectedAlgorithm]  = parameters{:};
%         Delta = 0.5;
%         N = 10;
%         receivPosit = [Delta*(-(N-1)/2:(N-1)/2)' zeros(N,1)];
%                 
%         % Transmitters
%         numSources = 1;
%         signalPowers = 1*ones(1, numSources);
%         
%         % Transmitters
%         transPosit = [dist.*sind(DoA), dist.*cosd(DoA)];
%         
%         processingEnvironment = struct;
%         switch selectedAlgorithm
%             case 'ULA'
%                 processingEnvironment.disable_ULA = false;
%             case 'SRA'
%                 processingEnvironment.sparse_ruler = [0,1,4,7,9];
%                 processingEnvironment.disable_ULA = true;
%             case 'NE'
%                 processingEnvironment.nested_elements = nestedArrayElements([2,3]);
%                 processingEnvironment.disable_ULA = true;
%             case 'CS'
%                 processingEnvironment.coprime_sampling = coprimeSamplingElements(2, 3);
%                 processingEnvironment.disable_ULA = true;
%         end
%                 processingEnvironment.investigated_theta = -90:90;
% 
%         
%         varargout = {DoA, transPosit, signalPowers, receivPosit, processingEnvironment, selectedAlgorithm};
%         
%     case '2sources'
%         [DoA1, DoA2, dist, selectedAlgorithm] = parameters{:};
%         Delta = 0.5;
%         N = 10;
%         receivPosit = [Delta*(-(N-1)/2:(N-1)/2)' zeros(N,1)];
%         
%         % Transmitters
%         numSources = 2;
%         signalPowers = 1*ones(1, numSources);
%         DoAs = [DoA1; DoA2];
%         transPosit = dist*[sind(DoAs), cosd(DoAs)];
%         
%         processingEnvironment = struct;
%         switch selectedAlgorithm
%             case 'ULA'
%                 processingEnvironment.disable_ULA = false;
%             case 'SRA'
%                 processingEnvironment.sparse_ruler = [0,1,4,7,9];
%                 processingEnvironment.disable_ULA = true;
%             case 'NE'
%                 processingEnvironment.nested_elements = nestedArrayElements([2,3]);
%                 processingEnvironment.disable_ULA = true;
%             case 'CS'
%                 processingEnvironment.coprime_sampling = coprimeSamplingElements(2, 3);
%                 processingEnvironment.disable_ULA = true;
%         end
%         processingEnvironment.investigated_theta = samplingVector( DoAs, -90, 90, 1, 10);
%         
%         varargout = {DoAs, transPosit, signalPowers, receivPosit, processingEnvironment, selectedAlgorithm};
%     
%     case 'DoAs_dist'
%         [DoAs, distances, selectedAlgorithm] = parameters{:};
%         DoAs = DoAs(:); distances = distances(:);
%         % DoAs and distances are column vectors with the same number of
%         % elements
%         Delta = 0.5;
%         N = 10;
%         receivPosit = [Delta*(-(N-1)/2:(N-1)/2)' zeros(N,1)];
%         
%         % Transmitters
%         numSources = numel(DoAs);
%         signalPowers = 1*ones(1, numSources);
%         transPosit = [distances.*sind(DoAs), distances.*cosd(DoAs)];
% 
%          processingEnvironment = struct;
%         switch selectedAlgorithm
%             case 'ULA'
%                 processingEnvironment.disable_ULA = false;
%             case 'SRA'
%                 processingEnvironment.sparse_ruler = [0,1,4,7,9];
%                 processingEnvironment.disable_ULA = true;
%             case 'NE'
%                 processingEnvironment.nested_elements = nestedArrayElements([2,3]);
%                 processingEnvironment.disable_ULA = true;
%             case 'CS'
%                 processingEnvironment.coprime_sampling = coprimeSamplingElements(2, 3);
%                 processingEnvironment.disable_ULA = true;
%         end
%         processingEnvironment.investigated_theta = samplingVector( DoAs, -90, 90, 1, 10);
%         
%         varargout = {DoAs, transPosit, signalPowers, receivPosit, processingEnvironment, selectedAlgorithm};
%     otherwise
%         error('formatInputs:notRecognizedOuptu', 'The option is not recognized')
% end
% end

