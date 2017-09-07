function varargout = getInfoMultipath(inputMethod, parameters, selectedOutput, varargin) % Varargin are extra output options
% getInfoMultipath( selectedSource, selectedOutput, selectedAlgorithm, maxDev, option, varargin )
% varargout = getInfoMultipath( mainDoA, signalPowers, multipath, elecDist, processingEnvironment, maxDev, selectedSource, selectedOutput, selectedAlgorithm )
% getInfoMultipath(inputMethod, parameters, outputOptions)

% Give form to inputs
[mainDoA, signalPowers, multipath, elecDist, numAntennas, processingEnvironment, selectedAlgorithm] = formatInputs(inputMethod, parameters);

% Calculate
[calcDoA_ULA, calcDoA_sp, calcDoA_ne, calcDoA_cs, ...
    calcHeight_ULA, calcHeight_sp, calcHeight_ne, calcHeight_cs, spectrumData] = calcDoAmultipath( mainDoA, signalPowers, multipath, elecDist, numAntennas, processingEnvironment );

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
    varargout{1} = getPeaksInfo(mainDoA, calcDoA, calcHeight, {selectedOutput}, varargin);
end

end

function varargout = formatInputs(option, parameters)
switch option
    case '1source1multipath'
        [mainDoA, secDoA, relAmplitude, phases, delay, selectedAlgorithm, numSourcesSearch, maxDev]  = parameters{:};
        multipath = {[secDoA; relAmplitude; phases; delay]};
        
        signalPowers = 1;
        elecDist = 0.5;
        
        % Calculate Investigated Theta
        numSources = numel(mainDoA);
        virtualDoAs = cell(1, numSources);
        for s = 1:numSources
            virtualDoAs{s} = [mainDoA(s), multipath{s}(1,:)];
        end
        virtualDoAs = cell2mat(virtualDoAs);

        
processingEnvironment = struct;
processingEnvironment.numSources= numSourcesSearch; % Number of sources that MUSIC is going to look for
switch selectedAlgorithm
    case 'ULA'
%         selectedAlgorithm = 'ULA';
processingEnvironment.disable_ULA = false;
    case 'SRA'
%         selectedAlgorithm = 'SparseRuler';
processingEnvironment.sparse_ruler = [0,1,4,7,9];
processingEnvironment.disable_ULA = true;
    case 'NE'
%         selectedAlgorithm = 'NestedElements';
processingEnvironment.nested_elements = nestedArrayElements([2,3]);
processingEnvironment.disable_ULA = true;
    case 'CS'
%         selectedAlgorithm = 'CoprimeSampling';
processingEnvironment.coprime_sampling = coprimeSamplingElements(2, 3);
processingEnvironment.disable_ULA = true;
end

% maxDev = 10;
minimum = max(-90, min(virtualDoAs) - maxDev);
        maximum = min(90, max(virtualDoAs) + maxDev);
        processingEnvironment.investigated_theta = samplingVector( virtualDoAs, minimum, maximum, 1, 10);

        varargout = {mainDoA, signalPowers, multipath, elecDist, 10, processingEnvironment, selectedAlgorithm};    
    otherwise    
        error('formatInputs:notRecognizedOuptu', 'The option is not recognized')
end
end

