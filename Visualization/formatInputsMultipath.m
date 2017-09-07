function varargout = formatInputsMultipath(option, parameters)
switch option
    case '1source1multipath'
        [mainDoA, secDoA, relAmplitude, phases, delay, numSourcesSearch, maxDev]  = parameters{:};
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
        processingEnvironment.disable_ULA = false;
        processingEnvironment.sparse_ruler = [0,1,4,7,9];
        processingEnvironment.nested_elements = nestedArrayElements([2,3]);
        processingEnvironment.coprime_sampling = coprimeSamplingElements(2, 3);
        processingEnvironment.investigated_theta = samplingVector( [mainDoA, secDoA], -90, 90, 1, 10);
        processingEnvironment.numSources= numSourcesSearch; % Number of sources that MUSIC is going to look for
        minimum = max(-90, min(virtualDoAs) - maxDev);
        maximum = min(90, max(virtualDoAs) + maxDev);
        processingEnvironment.investigated_theta = samplingVector( virtualDoAs, minimum, maximum, 1, 10);

        varargout = {mainDoA, signalPowers, multipath, elecDist, 10, processingEnvironment};    
    otherwise    
        error('formatInputs:notRecognizedOuptu', 'The option is not recognized')
end
end