function [varargout] = formatInputsNearField( inputOption, inputParameters, outputOption )
% [varargout] = formatInputsNearField( input_args )

processingEnvironment = struct;

switch inputOption
    case '1source'
        [DoA, dist, selectedAlgorithms]  = inputParameters{:};
        Delta = 0.5;
        N = 10;
        receivPosit = [Delta*(-(N-1)/2:(N-1)/2)' zeros(N,1)];
                
        % Transmitters
        numSources = 1;
        signalPowers = 1*ones(1, numSources);
        DoAs = DoA;
        transPosit = [dist.*sind(DoAs), dist.*cosd(DoAs)];

        processingEnvironment = processSelectedAlgorithm(selectedAlgorithms, processingEnvironment);
        
    case '2sources'
        [DoA1, DoA2, dist, selectedAlgorithms] = inputParameters{:};
        Delta = 0.5;
        N = 10;
        receivPosit = [Delta*(-(N-1)/2:(N-1)/2)' zeros(N,1)];
        
        % Transmitters
        numSources = 2;
        signalPowers = 1*ones(1, numSources);
        DoAs = [DoA1; DoA2];
        transPosit = dist*[sind(DoAs), cosd(DoAs)];

        processingEnvironment = processSelectedAlgorithm(selectedAlgorithms, processingEnvironment);
                
    case 'DoAs_dist'
        [DoAs, distances, selectedAlgorithms] = inputParameters{:};
        DoAs = DoAs(:); distances = distances(:);
        % DoAs and distances are column vectors with the same number of
        % elements
        Delta = 0.5;
        N = 10;
        receivPosit = [Delta*(-(N-1)/2:(N-1)/2)' zeros(N,1)];
        
        % Transmitters
        numSources = numel(DoAs);
        signalPowers = 1*ones(1, numSources);
        transPosit = [distances.*sind(DoAs), distances.*cosd(DoAs)];

        processingEnvironment = processSelectedAlgorithm(selectedAlgorithms, processingEnvironment);
        
    otherwise
        error('formatInputs:notRecognizedOuptu', 'The option is not recognized')
end

processingEnvironment.investigated_theta = samplingVector( DoAs, -90, 90, 1, 10);


switch outputOption
    case 'basic'
        varargout = {transPosit, signalPowers, receivPosit, processingEnvironment};
    case 'extended'
        varargout = {DoAs, transPosit, signalPowers, receivPosit, processingEnvironment, selectedAlgorithms};
    otherwise
        error('formatInputsNearField:outputOption', 'outputOption not recognized')
end

end

function processingEnvironment = processSelectedAlgorithm(selectedAlgorithms, processingEnvironment)

processingEnvironment.disable_ULA = true;
if isscalar(selectedAlgorithms)
    switch selectedAlgorithms
        case 'ULA'
            processingEnvironment.disable_ULA = false;
        case 'SRA'
            processingEnvironment.sparse_ruler = [0,1,4,7,9];
        case 'NE'
            processingEnvironment.nested_elements = nestedArrayElements([2,3]);
        case 'CS'
            processingEnvironment.coprime_sampling = coprimeSamplingElements(2, 3);
    end
else
    for k = 1:numel(selectedAlgorithms)
        switch selectedAlgorithms(k)
            case 'ULA'
                processingEnvironment.disable_ULA = false;
            case 'SRA'
                processingEnvironment.sparse_ruler = [0,1,4,7,9];
            case 'NE'
                processingEnvironment.nested_elements = nestedArrayElements([2,3]);
            case 'CS'
                processingEnvironment.coprime_sampling = coprimeSamplingElements(2, 3);
        end
    end
end
end

