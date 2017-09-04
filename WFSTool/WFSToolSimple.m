classdef WFSToolSimple < handle
    % One device, sinusoidal signals
    
    properties(SetAccess = public)
        changed % Variables changed since last update
        simplePerformance = true;
        
        % Noise source
        noiseSourceOrientation
        noiseSourceRadiationPattern
        virtual
        virtualVolume
        real
        realVolume
        noiseSourceChannelMapping
        signalsSpec % String that specifies the coefficient
        amplitude
        phase
        frequency
        
        % WFS Array
        WFSarrayRadiationPattern
        WFSarrayChannelMapping
        
        % Microphones
        receiverChannelMapping % Active channels of the receiver device

        % Loudspeakers
        loudspeakerChannelMapping
        
        % Calibration
        sourceCorr % Correction factor from adimensional variable to physic variable
        receiverCorr % Correction factor from physics magnitude of pressure to adimiensional units
                             
        % Experimental results
        pulseCoeffMat
        pulseLimits % Dimension 1 of pulseCoeffMat
        reprodChannelMapping % Dimension 2 of pulseCoeffMat
        reprodFrequencies % Dimension 3 of pulseCoeffMat
        recordedChannelMapping
        recordedSampleRate
        
        % Acoustic Paths
        % An acoustic path structure has 4 fields
        % - acousticPaths
        % - xChannelMapping
        % - yChannelMapping
        % - frequencies
        expAcPathStruct % Last experimental acoustic path calculated. If it's done well, it is supposed to be updated for the current physical scenario.
        loudspeakerAcPathStruct
        noiseSourceAcPathStruct
        WFSarrayAcPathStruct
        
        % Simulation results
        
        
        % Infrastructure variables
        fig
        player
        reprodPanel
        scenarioObj
        propPanel
        recordPanel
        simulObj
                
        % Other
        ax
        timeDisp
    end
    
    properties(Dependent)
        % Noise source
        noiseSourcePosition
   
        % WFS Array
        WFSarrayPosition % Dependent on scenarioObj
        WFSarrayOrientation % Dependent on scenarioObj
        
        % Microphones
        receiverPosition
        receiverOrientation
        receiverRadiationPattern
        
        % Loudspeakers
        loudspeakerPosition
        loudspeakerCoefficient
        loudspeakerOrientation
        loudspeakerRadiattionPattern
        loudspeakerFrequencies
        
        % Experiment
        recorded
        
        % Simulation
        simulField
        
        % Counting
        numNoiseSources
        numReceivers
        numSourcesWFSarray
        numLoudspeakers
    end
    
    properties(Constant)
        c = 340; % m/s
    end
    
    % Getters and setters
    methods
        % Counting
        function numNoiseSources = get.numNoiseSources(obj)
            % Number of noise sources
            numNoiseSources = numel(obj.signalsSpec);
        end   
        
        function numReceivers = get.numReceivers(obj)
            numReceivers = size(obj.receiverPosition, 1);
        end
        
        function numSourcesWFSarray = get.numSourcesWFSarray(obj)
            numSourcesWFSarray = obj.scenarioObj.numLoudspeakers;
        end
        
        function numLoudspeakers = get.numLoudspeakers(obj)
%             numLoudspeakers = obj.player.numChannels(1);
            numLoudspeakers = size(obj.loudspeakerCoefficient, 1);
        end
        
        % Noise source
        function noiseSourcePosition = get.noiseSourcePosition(obj)
            indActiveSour = find(obj.real | obj.virtual);
%             indReal = find(obj.real);
%             indRealInActive = ismember(indActiveSour, indReal);
%             
            noiseSourcePosition = zeros(obj.numNoiseSources, 3);
            noiseSourcePosition(indActiveSour, :) = obj.scenarioObj.sourcesPosition(indActiveSour, :);
        end
        
        function set.noiseSourcePosition(obj, value)
            numPos = size(value, 1);
            
            if numPos == obj.numNoiseSources
                obj.scenarioObj.setSourcePosition(value, (1:size(value, 1)));
            else
                warning('WFSToolSimple:WrongNumberNoiseSources', 'The number of positions must be equal to the number of sources');
            end
        end

        % WFS array
        function WFSarrayPosition = get.WFSarrayPosition(obj)
            WFSarrayPosition = obj.scenarioObj.loudspeakersPosition;
        end
        
        function set.WFSarrayPosition(obj, value)
            obj.scenarioObj.loudspeakersPosition = value;
        end

        function set.WFSarrayOrientation(obj, value)
            obj.scenarioObj.loudspeakersOrientation = simulator.rotVec2BroadsideVec(value);
        end

        function WFSarrayOrientation = get.WFSarrayOrientation(obj)
            WFSarrayOrientation = simulator.vec2rotVec(obj.scenarioObj.loudspeakersOrientation);
        end
        
        % Microphones
        function receiverPosition = get.receiverPosition(obj)
            receiverPosition = obj.scenarioObj.receiversPosition;
        end
        
        function set.receiverPosition(obj, value)
            obj.scenarioObj.setReceiverPosition(value, (1:size(value, 1)));
            obj.simulObj.recPositions = value;
        end
           
        function receiverOrientation = get.receiverOrientation(obj)
            receiverOrientation = obj.simulObj.recOrientations;
        end
        
        function set.receiverOrientation(obj, value)
            obj.simulObj.recOrientations = value;
        end
         
        function receiverRadiationPattern = get.receiverRadiationPattern(obj)
            receiverRadiationPattern = obj.simulObj.recRadPatFuns;
        end
        
        function set.receiverRadiationPattern(obj, value)
            obj.simulObj.recRadPatFuns = value;
        end        
        
        % Loudspeakers
        function loudspeakerPosition = get.loudspeakerPosition(obj)
            loudspeakerPosition = obj.simulObj.sourcePositions;
        end
        
        function set.loudspeakerPosition(obj, value)
            obj.simulObj.sourcePositions = value;
        end
        
        function loudspeakerCoefficient = get.loudspeakerCoefficient(obj)
            loudspeakerCoefficient = obj.simulObj.sourceCoefficients;
        end
        
        function set.loudspeakerCoefficient(obj, value)
            obj.simulObj.sourceCoefficients = value;
        end
        
        function loudspeakerOrientation = get.loudspeakerOrientation(obj)
            loudspeakerOrientation = obj.simulObj.sourceOrientations;
        end
        
        function set.loudspeakerOrientation(obj, value)
            obj.simulObj.sourceOrientations = value;
        end
        
        function loudspeakerRadiattionPattern = get.loudspeakerRadiattionPattern(obj)
            loudspeakerRadiattionPattern = obj.simulObj.radPatFuns;
        end
        
        function set.loudspeakerRadiattionPattern(obj, value)
            obj.simulObj.radPatFuns = value;
        end
        
        function loudspeakerFrequencies = get.loudspeakerFrequencies(obj)
            loudspeakerFrequencies = obj.simulObj.freq;
        end
        
        function set.loudspeakerFrequencies(obj, value)
            obj.simulObj.freq = value;
        end
       
        % Experimental
        function recorded = get.recorded(obj)
            recorded = obj.player.recorded{1};
        end
                
        % Simulation
        function simulField = get.simulField(obj)
            simulField = obj.simulObj.field;
        end
    end
    
    methods
        
        function obj = WFSToolSimple()
            fig = figure('Units', 'pixels', 'Position', [0 50 1200 600]);
            obj.fig = fig;
            
            obj.player = reproductorRecorder();
            obj.reprodPanel = reproductionPanel_noiseChannel(fig, [0.05, 0.6, 0.4, 0.4], @(action) obj.orderCallback(action));
            obj.scenarioObj = scenario(fig);
            obj.propPanel = propertiesPanel(fig, [0.05 0.1 0.4 0.2]);
            obj.recordPanel = recorderPanel(fig, [0.05 0.35 0.4 0.2]);
            obj.simulObj = simulator;
            obj.ax = obj.scenarioObj.ax;
            obj.simulObj.ax = obj.ax;
            colormap(obj.ax, 'gray')
            obj.ax.CLim = [-1 1];            
                        
            obj.changed = struct('virtual', false, 'real', false, 'signalsSpec', false);
            
            addlistener(obj.reprodPanel, 'updatedValues', @(~, evntData) obj.reprodPanelListener(evntData.type));
            addlistener(obj.recordPanel, 'updatedValues', @(~, evntData) obj.recordPanelListener(evntData.type));
            addlistener(obj.player, 'playingState', 'PostSet', @(~, eventData) obj.GUIenabling(eventData.AffectedObject.playingState));
            addlistener(obj.player, 'numChannels', 'PostSet', @(~, eventData) obj.changeScenario(eventData.AffectedObject.numChannels(1)));
            addlistener(obj.player, 'count', 'PostSet', @(~, eventData) obj.timeDisplay(eventData.AffectedObject.count*eventData.AffectedObject.frameDuration));
                                 
            obj.setNumNoiseSources(1);
            obj.setNumReceivers(1);
            obj.setNumWFSarraySources(2);
            
            obj.readFromReprodPanel();
            obj.changed.virtual = true;
            obj.changed.real = true;
            obj.changed.signalsSpec = true;
            obj.changed.activeReceivers = true;
            obj.changed.numNoiseSources = true;
            obj.updateEverything();          
            
            obj.simulObj.XnumPoints = 200;
            obj.simulObj.YnumPoints = 200;
            obj.simulateOnAxis();
            
            uicontrol(fig, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', [0.6 0.95 0.1 0.05], 'String', 'Simulate', 'Callback', @(hObject, eventData) obj.simulateOnAxis());
            uicontrol(fig, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', [0.75 0.95 0.1 0.05], 'String', 'Save', 'Callback', @(~, ~) obj.saveInformation());
            obj.timeDisp = uicontrol(fig, 'Style', 'text',...
                'Units', 'normalized', 'Position', [0.85 0.95 0.1 0.05], 'String', '');
        end
        
        function updateEverything(obj)
                       
            if obj.changed.numNoiseSources
                rightSize = [obj.numNoiseSources, 1];
                
                virtualRight = all(size(obj.reprodPanel.virtual) == rightSize);
                realRight = all(size(obj.reprodPanel.real) == rightSize);
                signalsRight = all(size(obj.signalsSpec) == rightSize);
                channelNumberRight = all(size(obj.noiseSourceChannelMapping) == rightSize);
                virtualVolumeRight = all(size(obj.virtualVolume) == rightSize);
                realVolumeRight = all(size(obj.virtualVolume) == rightSize);
                
                assert(virtualRight && realRight && signalsRight && channelNumberRight...
                    && virtualVolumeRight && realVolumeRight, 'WFSTool2:updateEverything', 'The signals specifications and the virtual and real flags must have the same size')
                
                comMat = WFSToolSimple.createCommutationMatrix(obj.virtual, obj.real);
                obj.player.setProps('comMatrix', comMat);
                
                obj.updateGUIConnectionsStuff_1Device();
                
                obj.updateSignalParameteres();
                obj.updateSignalProvidersVariables();
                obj.updateScenario(obj.real, obj.virtual, obj.reprodPanel.real, obj.reprodPanel.virtual);
                obj.updateDelayAndAttenFunctions();

                obj.changed.real = false;
                obj.changed.virtual = false;
                obj.changed.numNoiseSources = false;
                
            end
            
            if obj.changed.virtual || obj.changed.real                
                rightSize = [obj.numNoiseSources, 1];
                virtualRight = all(size(obj.virtual) == rightSize);
                realRight = all(size(obj.real) == rightSize);
                signalsRight = all(size(obj.signalsSpec) == rightSize);
                channelNumberRight = all(size(obj.noiseSourceChannelMapping) == rightSize);
                virtualVolumeRight = all(size(obj.virtualVolume) == rightSize);
                realVolumeRight = all(size(obj.virtualVolume) == rightSize);
                assert(virtualRight && realRight && signalsRight && channelNumberRight...
                    && virtualVolumeRight && realVolumeRight, 'WFSTool2:updateEverything', 'The signals specifications and the virtual and real flags must have the same size')
                
                obj.updateScenario(obj.real, obj.virtual, obj.reprodPanel.real, obj.reprodPanel.virtual);
                obj.updateDelayAndAttenFunctions();
               
                obj.changed.real = false;
                obj.changed.virtual = false;         
            end
            
            if obj.changed.signalsSpec
                obj.updateSignalParameteres();
                obj.updateSignalProvidersVariables();
            end
        end
        
        function setNumNoiseSources(obj, numNoiseSources)
            obj.virtual = true(numNoiseSources, 1);
            obj.real = true(numNoiseSources, 1);
            obj.signalsSpec = cell(numNoiseSources, 1);
            obj.noiseSourceChannelMapping = zeros(numNoiseSources, 1);
            obj.virtualVolume = ones(numNoiseSources, 1);
            obj.realVolume = ones(numNoiseSources, 1);
            obj.noiseSourceOrientation = repmat([1, 0, 0, 1], numNoiseSources, 1);
            obj.noiseSourceRadiationPattern = repmat({@ simulator.monopoleRadPat}, numNoiseSources, 1);
                        
            for k = 1:numNoiseSources
                obj.signalsSpec{k} = '';
            end
            
            obj.updateComMat();
            obj.updateSignalProvidersVariables();
            
            uniqueFreq = unique(obj.frequency);
            obj.noiseSourceAcPathStruct = struct('acousticPaths', zeros(obj.numReceivers, numNoiseSources, numel(uniqueFreq)), 'frequencies', uniqueFreq);
        end
        
        function setNumWFSarraySources(obj, numLoudspeakers)
                  
            s = WFSToolSimple.generateScenario(numLoudspeakers);
            obj.WFSarrayChannelMapping = (1:numLoudspeakers)';
            obj.WFSarrayRadiationPattern = repmat({@ simulator.monopoleRadPat}, numLoudspeakers, 1);
            uniqueFreq = unique(obj.frequency);
            obj.WFSarrayAcPathStruct = struct('acousticPaths', zeros(obj.numReceivers, numLoudspeakers, numel(uniqueFreq)), 'frequencies', uniqueFreq);
            
            % Reset noise sources and receivers
            sourcePosition = repmat([0 0 0], obj.numNoiseSources, 1);
            recPosition = repmat([0 0 0], obj.numReceivers, 1);
                        
            obj.scenarioObj.setScenario(sourcePosition, recPosition, s.loudspeakersPosition, s.loudspeakersOrientation, s.roomPosition);
                    
            obj.simulObj.XLim = [s.roomPosition(1), s.roomPosition(1) + s.roomPosition(3)];
            obj.simulObj.YLim = [s.roomPosition(2), s.roomPosition(2) + s.roomPosition(4)];
            delete(obj.simulObj.imag);
            
            obj.updateForcedDisabledLoudspeakers();
            
            obj.sourceCorr = ones(numLoudspeakers, 1);
        end
        
        function setNumReceivers(obj, numReceivers)
            obj.scenarioObj.setNumReceivers(numReceivers);
            obj.receiverPosition = obj.receiverPosition;
            obj.receiverOrientation = repmat([1, 0, 0, 1], numReceivers, 1);
            obj.receiverRadiationPattern = repmat({@ simulator.monopoleRadPat}, numReceivers);
            obj.receiverCorr = ones(numReceivers, 1);
            obj.receiverChannelMapping = (1:numReceivers)';
            uniqueFreq = unique(obj.frequency);
            obj.noiseSourceAcPathStruct = struct('acousticPaths', zeros(obj.numReceivers, obj.numNoiseSources, numel(uniqueFreq)), 'frequencies', uniqueFreq);
            obj.WFSarrayAcPathStruct = struct('acousticPaths', zeros(obj.numReceivers, obj.numSourcesWFSarray, numel(uniqueFreq)), 'frequencies', uniqueFreq);
        end
        
        function setVirtual(obj, value)
            obj.virtual = value;
            obj.changed.virtual = true;
        end
        
        function setReal(obj, value)
            obj.real = value;
            obj.changed.real = true;
        end
        
        function setSignalsSpec(obj, value)
            obj.signalsSpec = value;
            obj.changed.signalsSpec = true;
        end
        
        function setActiveReceivers(obj, value)
            obj.activeReceivers = value;
        end
        
        function calculateExperimentalAcousticPaths(obj)
            
            % Get variables           
            s = obj.exportInformation();
            sExp = s.Experiment;
                                    
            % Calculate the acoustic paths
            acPathStruct = struct;
            [acPathStruct.acousticPaths, acPathStruct.frequencies] = getAcousticPath( sExp.pulseLimits, sExp.frequencies, sExp.pulseCoefMat, sExp.recordedSignal, sExp.sampleRate);
            acPathStruct.xChannelMapping = sExp.channelMapping;
            acPathStruct.yChannelMapping = sExp.recordedChannelMapping;
            
            obj.expAcPathStruct = acPathStruct;
            
        end
                
        function WFScalculation(obj)
                    
            % We assume that ldspkrsCoef is a matrix of size
            % (numSourcesWFSarray x numNoiseSources). The (i, j)-th element
            % is the complex coefficient of the i-th WFS array element
            % destinated to synthesize the j-th noise source (with, of
            % course, the frequency corresponding to the j-th source). It
            % is not modified by the virtualVolume yet.
            
            % Unify the acoustic paths for the same frequencies
            uniqueFreq = unique(obj.frequency);
            WFSarrayAcPath = WFSToolSimple.tuneAcousticPaths(obj.WFSarrayAcPathStruct.acousticPaths, obj.WFSarrayAcPathStruct.frequencies, uniqueFreq);
            noiseSourcesAcPath = WFSToolSimple.tuneAcousticPaths(obj.noiseSourceAcPathStruct.acousticPaths, obj.noiseSourceAcPathStruct.frequencies, uniqueFreq);
            
            % Apply virtual flags and virtualVolume
            WFSarrayCoef = obj.getComplexCoeffWFS();
            WFSarrayCoef = WFSarrayCoef .* repmat(obj.virtualVolume', [obj.numSourcesWFSarray, 1]);
            WFSarrayCoef(:, ~obj.virtual) = 0;
            
            % Apply real flags and realVolume
            noiseSourceCoef = diag(obj.amplitude .* exp(1i*obj.phase) .* (-obj.realVolume));
            noiseSourceCoef(:, ~obj.real) = 0;   
            
            % Substitute data for the loudspeakers used as real noise
            % sources
            
                % Select the noise source values         
            nSrcChannelMapping_real = obj.noiseSourceChannelMapping(obj.real);
            nSrcCoef_real = noiseSourceCoef(obj.real, :);
            nSrcPos_real = obj.noiseSourcePosition(obj.real, :);
            nSrcOrient_real = obj.noiseSourceOrientation(obj.real, :); % simulator.vec2rotVec(repmat([0 0 1], [numReal, 1]));
            nSrcRadPat_real = obj.noiseSourceRadiationPattern(obj.real);
            nSrcAcPath_real = noiseSourcesAcPath(:, obj.real, :);
                             
            [ldspkrsChannelMapping, mapping] = combineIndices(obj.WFSarrayChannelMapping, nSrcChannelMapping_real);
            numLdsprks = numel(ldspkrsChannelMapping);
            
            ldspkrsCoef = zeros(numLdsprks, obj.numNoiseSources);
            ldspkrsPos = zeros(numLdsprks, 3);
            ldspkrsOrient = zeros(numLdsprks, 4);
            ldspkrsRadPat = cell(numLdsprks, 1);
            ldspkrsAcPath = zeros(obj.numReceivers, numLdsprks, numel(uniqueFreq));
            
            ldspkrsCoef(mapping(1).destinationInd, :) = WFSarrayCoef(mapping(1).originInd, :);
            ldspkrsPos(mapping(1).destinationInd, :) = obj.WFSarrayPosition(mapping(1).originInd, :);
            ldspkrsOrient(mapping(1).destinationInd, :) = obj.WFSarrayOrientation(mapping(1).originInd, :);
            ldspkrsRadPat(mapping(1).destinationInd) = obj.WFSarrayRadiationPattern(mapping(1).originInd);
            ldspkrsAcPath(:, mapping(1).destinationInd, :) = WFSarrayAcPath(:, mapping(1).originInd, :);
            
            ldspkrsCoef(mapping(2).destinationInd, :) = nSrcCoef_real(mapping(2).originInd, :);
            ldspkrsPos(mapping(2).destinationInd, :) = nSrcPos_real(mapping(2).originInd, :);
            ldspkrsOrient(mapping(2).destinationInd, :) = nSrcOrient_real(mapping(2).originInd, :);
            ldspkrsRadPat(mapping(2).destinationInd) = nSrcRadPat_real(mapping(2).originInd);
            ldspkrsAcPath(:, mapping(2).destinationInd, :) = nSrcAcPath_real(:, mapping(2).originInd, :);
                       
            % Prepare acoustic paths structure
            ldspkrsAcPathStruct = struct('acousticPaths', ldspkrsAcPath, 'frequencies', uniqueFreq); % completar
            
            % Set the variables in the simulation object
            obj.loudspeakerPosition = ldspkrsPos;
            obj.loudspeakerCoefficient = ldspkrsCoef;
            obj.loudspeakerOrientation = ldspkrsOrient;
            obj.loudspeakerRadiattionPattern = ldspkrsRadPat; % repmat({@(x) simulator.monopoleRadPat(x)}, [obj.numLoudspeakers, 1]);
            obj.loudspeakerFrequencies = obj.frequency;
            obj.loudspeakerChannelMapping = ldspkrsChannelMapping;
            obj.loudspeakerAcPathStruct = ldspkrsAcPathStruct;

        end
              
        function reproduceSignalFunction(obj, signalFunction, signalSampleRate, channelMapping)
               
            % Save ReproductorRecorder parameters
            real_comMat = obj.player.comMatrix;
            writingDrivers = obj.player.driver;
            writingDevices = obj.player.device;
            
            % Change state
            obj.player.setProps('comMatrix', true);
            obj.player.setProps('mode', originType('func'), 1);
            obj.player.setProps('FsGenerator', signalSampleRate, 1);
            obj.player.setProps('enableProc', false);
            obj.player.setProps('driver', writingDrivers{1}, 1);
            obj.player.setProps('device', writingDevices{1}, 1);
            obj.player.setProps('signalFunc', signalFunction, 1);
            if nargin == 4
                obj.player.setProps('defaultChannelMapping', false, 1);
                obj.player.setProps('channelMapping_player', channelMapping, 1);
            else
                channelMapping = [];
            end

            obj.player.executeOrder('play');
                        
            % Save information about reproduction
            obj.reprodChannelMapping = channelMapping;
            obj.recordedChannelMapping = obj.player.channelMapping_recorder{1};
            obj.recordedSampleRate = obj.player.Fs_recorder(1);
            
            % Turn everything to the previous state
            obj.player.setProps('enableProc', true);
            obj.player.setProps('comMatrix', real_comMat);
            obj.updateSignalProvidersVariables();
            for k = 1:obj.player.numPlayers
                obj.player.setProps('driver', writingDrivers{k}, k);
                obj.player.setProps('device', writingDevices{k}, k);
            end        
            
        end
        
        function reproduceAndRecord(obj, option, varargin)
            SampleRate = 44100;
            frequencies = obj.loudspeakerFrequencies;
            
            [ pulseCoefMat, pulseLim ] = coefficients2pulseSignalParameters( obj.loudspeakerCoefficient, frequencies, SampleRate, option, varargin{:} );
            signalFunc = @(startSample, endSample) pulseCoefMat2signal(frequencies, pulseCoefMat, pulseLim, SampleRate, startSample, endSample, 'sample');           
%             signalFunc = @(startSample, endSample) coefficients2signal( sourcesCoef, frequencies, SampleRate, startSample, endSample, false, option, varargin{:});

            obj.reproduceSignalFunction(signalFunc, SampleRate, obj.loudspeakerChannelMapping);
                           
            obj.pulseCoeffMat = pulseCoefMat;
            obj.pulseLimits = pulseLim/SampleRate;
            obj.reprodFrequencies = frequencies;
            
        end        
        
        function simulate(obj)
            % Set the acoustic path
            simulAcPath = WFSToolSimple.tuneAcousticPaths(obj.loudspeakerAcPathStruct.acousticPaths, obj.loudspeakerAcPathStruct.frequencies, obj.loudspeakerFrequencies);
            obj.simulObj.acPath = simulAcPath;

            % Simulate
            obj.simulObj.updateField();
            
            % Draw
            obj.simulObj.drawingOption = 'scatter';
            obj.simulObj.draw();
            
        end        
        
        function reproduceAndRecordForAcousticPaths(obj)
            SampleRate = 44100;
            numChannels = obj.numLoudspeakers;
            frequencies = obj.frequency;
            numFreq = numel(frequencies);
            
            % Set coefficients for calibration
            calCoeff = 0.4*ones(numChannels, numFreq);
            
            % Reproduce
            [ pulseCoefMat, pulseLim ] = coefficients2pulseSignalParameters( calCoeff, frequencies, SampleRate, 'preludeAndMain', 1, 1, 2 );
            signalFunc = @(startSample, endSample) pulseCoefMat2signal(pulseCoefMat, pulseLim, frequencies, SampleRate, startSample, endSample, 'type_marker', 'sample');          
            obj.reproduceSignalFunction(signalFunc, SampleRate);
            
            % Save information about the reproduced signal
            obj.pulseCoeffMat = pulseCoefMat;
            obj.pulseLimits = pulseLim/SampleRate;
            obj.reprodFrequencies = frequencies;
        end
        
        function calibrate(obj)
                        
            % Calculate experimental acoustic paths
            obj.reproduceAndRecordForAcousticPaths();
            obj.calculateExperimentalAcousticPaths();
            
            % Simulate acoustic path
            obj.calculateSimulatedAcousticPaths();
            
            % Get the ratio between them
            rat = obj.expAcPath./obj.simulAcPath;
            
            % Set source coefficient corrections (Best rank 1 approximation)
            [U, S, V] = svd(rat);
            obj.sourceCorr = U(:, 1);
            obj.receiverCorr = V(:, 1)*S(1);
            
        end
        
        function position = calculatePosition(obj)
            
            % Create pulse coefficient matrix that uses different
            % frequencies for each loudspeaker in order to calculate the
            % distance
            SampleRate = 44100;
            freqVec = [600, 610]';
            numFreq = numel(freqVec);
            numChannels = obj.numLoudspeakers;
            coefMat = 0.1*ones(numChannels, numFreq);         
            [ pulseCoefMat, pulseLim ] = coefficients2pulseSignalParameters( coefMat, freqVec, SampleRate, 'prelude', 0.25, 0.25, 1 );
            M = zeros(numFreq, numChannels);
            for f = 1:numFreq
                M(f, :) = (1:numChannels) + (f-1)*numChannels;
            end
            pulseCoefMat = pulseCoefMat(M(:), :, :);
            
            % Reproduce and Record
            signalFunc = @(startSample, endSample) pulseCoefMat2signal(freqVec, pulseCoefMat, pulseLim, SampleRate, startSample, endSample, 'sample');          
            obj.reproduceSignalFunction(signalFunc, SampleRate);
            
            % Save information about the reproduced signal
            obj.pulseCoeffMat = pulseCoefMat;
            obj.pulseLimits = pulseLim/SampleRate;
            obj.reprodFrequencies = freqVec;
            
            % Calculate the coefficients of the received signal
            s = obj.exportInformation();
            sExp = s.Experiment;
            acPath = getAcousticPath( freqVec, sExp.pulseCoefMat, sExp.pulseLimits, sExp.recordedSignal, sExp.recordedSignal_SampleRate);
            
            % Calculate distances
            phases = -angle(acPath);
            deltaPhase = diff(phases, 1, 3);
            deltaPhase = mod(deltaPhase, 2*pi);
            deltaFreq = diff(freqVec);
            deltaFreqMat = repmat(permute(deltaFreq, [2, 3, 1]), [numChannels, obj.numReceivers, 1]);
            dist = deltaPhase./deltaFreqMat*obj.c/(2*pi);
            dist = mean(dist, 3);
            
            % Calculate position
            refPoints = s.sScen.sourcesPos;
            position = zeros(obj.numReceivers, 3);
            for rec = 1:obj.numReceivers
                position(rec, :) = multilateration(refPoints, dist);
            end
            
        end
        
        function updateReprodPanelBasedOnVariables(obj)
            obj.reprodPanel.setNumSignals(obj.numNoiseSources);
            obj.reprodPanel.setRealFlags(obj.real);
            obj.reprodPanel.setVirtualFlags(obj.virtual);
            obj.reprodPanel.setChannelNumber(obj.noiseSourceChannelMapping);
            obj.updateSignalSpecFromParameters();
            obj.reprodPanel.setSignals(obj.signalsSpec);
        end
        
    end
    
    % Export and import methods
    methods
               
        function saveInformation(obj)
            defaultName = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
            [FileName,PathName, ~] = uiputfile('*.mat', 'Save Information', ['../Data/', defaultName]);
            
            if FileName ~= 0
                
                [~, ~, ext] = fileparts(FileName);
                if isempty(ext)
                    FileName = [FileName, '.mat'];
                end
                
                s = obj.exportInformation();
                save([PathName, FileName], 's');            
            end
        end
               
        function s = exportInformation(obj)
                        
            % Noise source variables
            s.NoiseSource = obj.getNoiseSourceVariables();
            
            % WFS array variables
            s.WFSarray = obj.getWFSarrayVariables();
            
            % Receiver variables
            s.Receiver = obj.getReceiverVariables();
            
            % Loudspeaker variables
            s.Loudspeaker = obj.getLoudspeakerVariables();
            
            % Simulation Results Variables
            s.Simulation = obj.getSimulationResultVariables();
            
            % Experimental Results Variables
            s.Experiment = obj.getExperimentalResultVariables();

        end
        
        function s = getExperimentalResultVariables(obj)
                        
            s.recordedSignal = obj.recorded;
            s.sampleRate = obj.recordedSampleRate;
            s.recordedChannelMapping = obj.recordedChannelMapping;
            s.pulseCoefMat = obj.pulseCoeffMat;
            s.pulseLimits = obj.pulseLimits;
            s.channelMapping = obj.reprodChannelMapping;
            s.frequencies = obj.reprodFrequencies;
            s.acPathStruct = obj.expAcPathStruct; % Sometimes it won't have been calculated
                   
        end
        
        function s = getSimulationResultVariables(obj)
            s.simulatedField = obj.simulField;
        end
        
        function s = getNoiseSourceVariables(obj)
            % Noise source
            s.noiseSourceOrientation = obj.noiseSourceOrientation;
            s.noiseSourceRadiationPattern = obj.noiseSourceRadiationPattern;
            s.virtual = obj.virtual;
            s.virtualVolume = obj.virtualVolume;
            s.real = obj.real;
            s.realVolume = obj.realVolume;
            s.noiseSourceChannelMapping = obj.noiseSourceChannelMapping;
            s.amplitude = obj.amplitude;
            s.phase = obj.phase;
            s.frequency = obj.frequency;
            s.noiseSourcePosition = obj.noiseSourcePosition;
        end
        
        function s = getWFSarrayVariables(obj)
            s.WFSarrayRadiationPattern = obj.WFSarrayRadiationPattern;
            s.WFSarrayChannelMapping = obj.WFSarrayChannelMapping;
            s.WFSarrayPosition = obj.WFSarrayPosition;
            s.WFSarrayOrientation = obj.WFSarrayOrientation;
        end
        
        function s = getReceiverVariables(obj)
            s.receiverOrientation = obj.receiverOrientation;
            s.receiverRadiationPattern = obj.receiverRadiationPattern;
            s.receiverChannelMapping = obj.receiverChannelMapping;
            s.receiverPosition = obj.receiverPosition;
        end
        
        function s = getLoudspeakerVariables(obj)
            s.loudspeakerChannelMapping = obj.loudspeakerChannelMapping;
            s.loudspeakerPosition = obj.loudspeakerPosition;
            s.loudspeakerCoefficient = obj.loudspeakerCoefficient;
            s.loudspeakerOrientation = obj.loudspeakerOrientation;
            s.loudspeakerRadiattionPattern = obj.loudspeakerRadiattionPattern;
            s.loudspeakerFrequencies = obj.loudspeakerFrequencies;
        end     
        
    end
    
    % WFS actions
    methods
        
        function WFS2realRatio(obj)
            
            realChannels = obj.noiseSourceChannelMapping(obj.real);
            realFlag = ismember(obj.loudspeakerChannelMapping, realChannels);
            WFSflag = ~realFlag;
            
            xLim = [1, 2]; yLim = [2, 4];
            x = linspace(xLim(1), xLim(2), 100);
            y = linspace(yLim(1), yLim(2), 100);
            z = 0;
            [X, Y, Z] = ndgrid(x, y, z);
            testPoints = [X(:), Y(:), Z(:)];
            
            U = mean(obj.simulObj.calculate(testPoints), 2); % Mean along test points
            U_WFS = sum(U(WFSflag, :, :), 1);
            U_real = sum(U(realFlag, :, :), 1);
            
            correction = -U_real./U_WFS; % (1 x numTestPoints x numNoiseSources)
            
            correction = permute(mean(correction, 2), [1 3 2]);
            
            obj.loudspeakerCoefficient(WFSflag, :) = obj.loudspeakerCoefficient(WFSflag, :).*repmat(correction, [sum(WFSflag), 1]);
                        
        end
             
        function nullField_(obj)
                                  
            coefficients_new = WFSToolSimple.nullField(obj.loudspeakerCoefficient,...
                obj.loudspeakerChannelMapping, obj.loudspeakerFrequencies,...
                obj.noiseSourceChannelMapping(obj.real), obj.expAcPathStruct);
            
            obj.loudspeakerCoefficient = coefficients_new;
        end
                 
    end
    
    
    methods(Access = private)
        
        function timeDisplay(obj, time)
            obj.timeDisp.String = num2str(time);
        end
        
        function setDefaultValues(obj)
            obj.setNumNoiseSources(1);
            obj.setNumReceivers(1);
            
        end
        
        function readFromReprodPanel(obj)
            obj.signalsSpec = obj.reprodPanel.signals;
            obj.virtual = obj.reprodPanel.virtual;
            obj.real = obj.reprodPanel.real;
            obj.noiseSourceChannelMapping = obj.reprodPanel.channelNumber;
            obj.virtualVolume = obj.reprodPanel.virtualVolume;
            obj.realVolume = obj.reprodPanel.realVolume;
        end
        
        function GUIenabling(obj, newPlayingState)
            if newPlayingState == playingStateClass('stopped')
%                 obj.reprodPanel.enableGUI();
                obj.propPanel.enableGUI();
            else
%                 obj.reprodPanel.disableGUI();
                obj.propPanel.disableGUI();
            end
        end
        
        function reprodPanelListener(obj, name)
            switch name
                case 'signals'
                    obj.signalsSpec = obj.reprodPanel.signals;
                    obj.changed.signalsSpec = true;
                case 'numSources'
                    obj.signalsSpec = obj.reprodPanel.signals;
                    obj.noiseSourceChannelMapping = obj.reprodPanel.channelNumber;
                    obj.virtualVolume = obj.reprodPanel.virtualVolume;
                    obj.realVolume = obj.reprodPanel.realVolume;
                    
                    obj.changed.numNoiseSources = true;
                case 'virtual'
                    obj.changed.virtual = true;
                case 'real'
                    obj.changed.real = true;
                case 'channelNumber'
                    obj.noiseSourceChannelMapping = obj.reprodPanel.channelNumber;
                    obj.updateForcedDisabledLoudspeakers();
                case 'virtualVolume'
                    obj.virtualVolume = obj.reprodPanel.virtualVolume;
                case 'realVolume'
                    obj.realVolume = obj.reprodPanel.realVolume;
            end
            
            obj.updateEverything();
        end
        
        function recordPanelListener(obj, ~)
            numRec = numel(obj.recordPanel.activeChannels);
            obj.setNumReceivers(numRec);
            obj.receiverChannelMapping = obj.recordPanel.activeChannels;
            obj.player.setProps('channelMapping_recorder', obj.receiverChannelMapping, 1);
        end
        
        function updateScenario(obj, oldReal, oldVirtual, newReal, newVirtual)
            oldN = numel(oldReal);
            newN = numel(newReal);
            if oldN ~= newN
                obj.scenarioObj.setNumSources(newN);
            else
                indActive = find(oldReal | oldVirtual);
                indNewActive = find(newReal | newVirtual);
                
                [a, b] = ismember(indNewActive, indActive);
                newPositions = zeros(numel(indNewActive), 3);
                newPositions(a, :) = obj.scenarioObj.sourcesPosition(b(a), :);
                
                obj.scenarioObj.setNumSources(size(newPositions, 1));
                obj.scenarioObj.setSourcePosition(newPositions, 1:numel(indNewActive));
            end
            
            obj.real = newReal;
            obj.virtual = newVirtual;
            
            obj.updateForcedDisabledLoudspeakers();
                       
        end
        
        function updateComMat(obj)            
            % Update commutation matrix
            comMat = WFSToolSimple.createCommutationMatrix(obj.virtual, obj.real);
            obj.player.setProps('comMatrix', comMat);
            
            obj.updateSignalParameteres();
            obj.updateSignalProvidersVariables();
            obj.updateProcessorVariables_noiseChannel();
            obj.updateGUIConnectionsStuff();
            
            obj.changed.real = false;
            obj.changed.virtual = false;
            obj.changed.numNoiseSources = false;
        end

        function updateSignalParameteres(obj)
            
            obj.amplitude = zeros(obj.numNoiseSources, 1);
            obj.phase = zeros(obj.numNoiseSources, 1);
            obj.frequency = zeros(obj.numNoiseSources, 1);
            
            for k = 1:obj.numNoiseSources
                signalSpec = obj.signalsSpec{k};
                
                % Is a complex number?
                param = regexp(signalSpec, 'A:(?<Amplitude>(\d+\.\d+|\d+)) Ph:(?<Phase>(\d+\.\d+|\d+)) f:(?<Frequency>(\d+\.\d+|\d+))', 'names');
                if isempty(param)
                    % It is not a complex number
                    % As default, treat it as a tone
                    obj.amplitude(k) = 0;
                    obj.phase(k) = 0;
                    obj.frequency(k) = 1;
                else
                    % It is a complex number, set the
                    % properties
                    obj.amplitude(k) = str2double(param.Amplitude);
                    obj.phase(k) = str2double(param.Phase);
                    obj.frequency(k) = str2double(param.Frequency);
                end
            end
        end
        
        function updateSignalSpecFromParameters(obj)
            signSpec = cell(obj.numNoiseSources, 1);
            for k = 1:obj.numNoiseSources
                signSpec{k} = sprintf('A:%.2f Ph:%.2f f:%.1f', obj.amplitude(k), obj.phase(k), obj.frequency(k));
            end
            obj.signalsSpec = signSpec;
        end
        
        function updateSignalProvidersVariables(obj)
            % Signal providers need to have some parameters specified based
            % on the signals specifications            
            % Assign signals
            for k = 1:obj.numNoiseSources
                obj.player.setProps('mode', originType('sinusoidal'), k);
                obj.player.setProps('amplitude', obj.amplitude(k), k);
                obj.player.setProps('phase', obj.phase(k), k);
                obj.player.setProps('frequency', obj.frequency(k), k);
            end
            
            obj.changed.signalsSpec = false;
        end
                
        function updateDelayAndAttenFunctions(obj)
            % What processors exist?
            indExist = find(obj.player.comMatrix);
            
            % Assign delay and attenuation functions            
            for k = 1:numel(indExist)
                index = indExist(k);
                
                delayFun = @() obj.delayFunction(index);
                obj.player.setProps('getDelayFun', delayFun, [index, 1]);
                
                attenFun = @() obj.attenuationFunction(index);
                obj.player.setProps('getAttenFun', attenFun, [index, 1]);
            end
        end
        
        function updateProcessorVariables_noiseChannel(obj)
            % The processors need the getDelayFun and getAttenFun
            % functions, and these are connected to the scenario
            
            % Update scenario based on the virtual sources
            activeSources = obj.virtual | obj.real;
            obj.scenarioObj.setNumSources(sum(activeSources));
            obj.updateForcedDisabledLoudspeakers();
            
            obj.updateDelayAndAttenFunctions();
            
        end

        function updateForcedDisabledLoudspeakers(obj)
            disabledLoudspeakers = false(obj.numSourcesWFSarray, 1);
            % Assign delay and attenuation functions
            for k = 1:obj.numNoiseSources
                if obj.real(k)
                    chann = obj.noiseSourceChannelMapping(k);
                    if chann > 0 && chann <= obj.numSourcesWFSarray
                        disabledLoudspeakers(chann) = true;
                    end
                end
            end
            obj.scenarioObj.setForcedDisabledLoudspeakers(disabledLoudspeakers);
        end
                
        function updateGUIConnectionsStuff(obj)
            % The devices need to have some variables specified:
            % driver and device
            % It creates a panel for configuring this variables as well as
            % the frame duration and the volume of each writing device, and
            % restart this properties to some defaults
            
            % Update player GUI controls based on the real sources
            %             indReal = [1, find(any(obj.player.comMatrix(:, 2:end), 1))+1];
            indWriters = find(any(obj.player.comMatrix, 1));
            obj.updateReproductorRecorderGUI_Writers(indWriters);
            
            indRec = 1:obj.player.numRecorders;
            obj.updateReproductorRecorderGUI_Readers(indRec);
        end    
        
        function updateGUIConnectionsStuff_1Device(obj)
            % The devices need to have some variables specified:
            % driver and device
            % It creates a panel for configuring this variables as well as
            % the frame duration and the volume of each writing device, and
            % restart this properties to some defaults
            
            % Update player GUI controls based on the real sources
            obj.updateReproductorRecorderGUI_Writers(1); % Only for WFSToolSimple. One device always active: WFS Array.
            
            obj.updateReproductorRecorderGUI_Readers(1);
            
        end
        
        function updateReproductorRecorderGUI_Writers(obj, indWritingDevices)
            setFrameDurationFunc = @(frameDuration) obj.player.setProps('frameDuration', frameDuration);
            setVolumeFunc = @(volume, index) obj.setVolume(volume, indWritingDevices(index));
            setDeviceFunc = @(device, index) obj.player.setProps('device', device, indWritingDevices(index));
            setDriverFunc = @(driver, index) obj.player.setProps('driver', driver, indWritingDevices(index));
            getAvailableDevicesFunc = @(index) obj.player.player{indWritingDevices(index)}.getAudioDevices();
            labels = cell(numel(indWritingDevices), 1);
            for k = 1:numel(indWritingDevices)
                labels{k} = sprintf('Device %d', indWritingDevices(k)-1);
            end
            if ~isempty(indWritingDevices) && indWritingDevices(1) == 1
                labels{1} = 'WFS Array';
            end
            
            obj.propPanel.setFunctionsWriters(setFrameDurationFunc, setVolumeFunc, setDeviceFunc, setDriverFunc, getAvailableDevicesFunc, labels);
            
        end
        
        function updateReproductorRecorderGUI_Readers(obj, indReadingDevices)
            setDeviceFunc = @(device, index) obj.player.setProps('device_recorder', device, indReadingDevices(index));
            setDriverFunc = @(driver, index) obj.player.setProps('driver_recorder', driver, indReadingDevices(index));
            getAvailableDevicesFunc = @(index) obj.player.recorder{indReadingDevices(index)}.getAudioDevices();
            labels = cell(numel(indReadingDevices), 1);
            for k = 1:numel(indReadingDevices)
                labels{k} = sprintf('Recorder %d', indReadingDevices(k));
            end
                        
            obj.propPanel.setFunctionsReaders(setDeviceFunc, setDriverFunc, getAvailableDevicesFunc, labels);
            
        end
      
        function delays = delayFunction(obj, index)
                      
            delays = obj.getDelays(index);
            delayShift = angle(obj.sourceCorr)/(2*pi*obj.frequency(index));
            delays = delays + delayShift;
            
        end
        
        function attenuations = attenuationFunction(obj, index)
            
            attenuations = obj.getAttenuations(index)./abs(obj.sourceCorr);
            
        end
             
        function delays = getDelays(obj, index)
            % Find index of source in scenarioObj. Only the sources that
            % are active (virtual, real or both flags set to true) are in
            % scenarioObj            
            activeSources = find(obj.virtual | obj.real);
            ind = (activeSources == activeSources(index));
            
            % Virtual
            if obj.virtual(index)
                delays = obj.scenarioObj.delays(:, ind);
            else
                delays= zeros(obj.numSourcesWFSarray, 1);
            end
            
            % If the real flag is set to true, the real value overrides the
            % virtual value for the selected loudspeaker
            if obj.real(index)              
                delays(obj.noiseSourceChannelMapping(index)) = 0;
            end
        end
        
        function delays = getDelaysWFS(obj, index)
                        
            if obj.virtual(index)
                % Find index of source in scenarioObj. Only the sources that
                % are active (virtual, real or both flags set to true) are in
                % scenarioObj
                activeSources = find(obj.virtual | obj.real);
                ind = (activeSources == activeSources(index));
                
                delays = obj.scenarioObj.delays(:, ind);
            else
                delays= zeros(obj.numSourcesWFSarray, 1);
            end
            
        end
        
        function attenuations = getAttenuations(obj, index)
            % Find index of source in scenarioObj. Only the sources that
            % are active (virtual, real or both flags set to true) are in
            % scenarioObj            
            activeSources = find(obj.virtual | obj.real);
            ind = (activeSources == activeSources(index));
            
            % Virtual
            if obj.virtual(index)
                attenuations = obj.virtualVolume(index)*obj.scenarioObj.attenuations(:, ind);
            else
                attenuations = zeros(obj.numSourcesWFSarray, 1);
            end
            
            % If the real flag is set to true, the real value overrides the
            % virtual value for the selected loudspeaker
            if obj.real(index) 
                attenuations(obj.noiseSourceChannelMapping(index)) = -obj.realVolume(index);
            end
        end
        
        function attenuations = getAttenuationsWFS(obj, index)
                    
            if obj.virtual(index)
                % Find index of source in scenarioObj. Only the sources that
                % are active (virtual, real or both flags set to true) are in
                % scenarioObj
                activeSources = find(obj.virtual | obj.real);
                ind = (activeSources == activeSources(index));
                
                attenuations = obj.scenarioObj.attenuations(:, ind);
            else
                attenuations = zeros(obj.numSourcesWFSarray, 1);
            end
            
        end
        
        function complexCoeff = getComplexCoeffWFS(obj)
            complexCoeff = zeros(obj.numSourcesWFSarray, obj.numNoiseSources);
            for k = 1:obj.numNoiseSources
                delays = obj.getDelaysWFS(k);
                atten = obj.getAttenuationsWFS(k);
                                
                amp = obj.amplitude(k);
                pha = obj.phase(k);
                sourceCoef = amp * exp(1i*pha);
                
                complexCoeff(:, k) = sourceCoef * atten .* exp(-1i*2*pi*obj.frequency(k)*delays);
            end
        end
               
        function complexCoeff = getComplexCoeff(obj)
            numFreq = obj.numNoiseSources;
            complexCoeff = zeros(obj.numLoudspeakers, numFreq);
            for k = 1:numFreq
                delays = obj.getDelays(k);
                atten = obj.getAttenuations(k);
                
                amp = obj.amplitude(k);
                pha = obj.phase(k);
                sourceCoef = amp * exp(1i*pha);
                
                complexCoeff(:, k) = sourceCoef * atten .* exp(-1i*2*pi*obj.frequency(k)*delays);
            end
        end
        
        function setVolume(obj, volume, indPlayer)
            comMatCoef = obj.player.comMatrixCoef;
            comMatCoef(:, indPlayer) = volume;
            obj.player.setProps('comMatrixCoef', comMatCoef);
        end
              
        function orderCallback(obj, order)
            % Based on the user order and the state of the player, a
            % command for the player is created
            state = obj.player.playingState;
            
            switch order
                case 'play'
                    % Start track again.
                    obj.player.executeOrder('stop');
                    if obj.simplePerformance
%                         obj.WFScalculation();
                        obj.reproduceAndRecord();
%                         obj.simulate();
                    else
                        obj.player.executeOrder('play');
                    end
                    
                    
                case 'stop'
                    % Stop
                    obj.player.executeOrder('stop');
                case 'pause'
                    % Pause
                    switch state
                        case playingStateClass('playing')
                            obj.player.executeOrder('pause');
                        case playingStateClass('stopped')
                            % Do nothing
                        case playingStateClass('paused')
                            obj.player.executeOrder('resume');
                    end
            end
            
        end
                  
        function simulateOnAxis(obj)
            % Configure simulator object
            obj.WFScalculation();
            
            % Simulate
            obj.simulObj.generateMeasurePoints();
            obj.simulObj.updateTheoricAcousticPathsImage();
            obj.simulObj.updateFieldImage();
            
            % Draw
            obj.simulObj.drawingOption = 'image';
            obj.simulObj.draw();
        end
        
    end
    
    methods(Static)
        
        function comMat = createCommutationMatrix(virtual, ~)
%             comMat = virtual | real;
            comMat = true(numel(virtual), 1);
        end
        
        function s = generateScenario(numLoudspeakers)
            if ~ismember(numLoudspeakers, [0 2 96])
                numLoudspeakers = 0;
            end
            
            switch numLoudspeakers
                case 0

                    loudspeakersPosition = double.empty(0,3);
                    loudspeakersOrientation = double.empty(0,3);
                    roomPosition = [0 0 1 1];   
                    
                case 2
                    loudspeakersPosition = [-0.1 0 0; 0.1 0 0];
                    loudspeakersOrientation = [1 0 0; -1 0 0];
                    roomPosition = [-2, -2, 4, 4];
                   
                case 96
                    d = 0.18; % Separation between two contiguous loudspeakers. Size of one loudspeaker
                    nb = 8; % Bottom and upper sides of the octogon (2 sides)
                    nd = 8; % Diagonal sides of the octogon (4 sides)
                    nl = 24; % Lateral side of the octogon (2 sides)
                    betabd = 45; % Deviation angle between bottom/upper and diagonal sides
                    
                    [ x, y, alfa ] = octogon(d, nb, nd, nl, betabd);
                    z = zeros(numel(x), 1);
                    loudspeakersPosition = [x, y, z];
                    loudspeakersOrientation = [cosd(alfa), sind(alfa), zeros(numel(alfa), 1)];
                    
                    xmin = min(x); xmax = max(x); ymin = min(y); ymax = max(y);
                    xDim = xmax - xmin; yDim = ymax - ymin;
                    xmargin = 0.2 * xDim; ymargin = 0.2 * yDim;
                    roomPosition = [xmin - xmargin, ymin - ymargin, xDim + 2*xmargin, yDim + 2*ymargin];
                    
                otherwise
                    warning('Wrong number of output channels. There is not possible scenario for that case')
            end
            
            s.loudspeakersPosition = loudspeakersPosition;
            s.loudspeakersOrientation = loudspeakersOrientation;
            s.roomPosition = roomPosition;
            
        end
        
        function acPath_query = tuneAcousticPaths(acPath, acPathFrequencies, queryFrequencies)
            % Based on the acoustic paths measured (or simulated) at some
            % specific frequencies (acPathFrequencies), get the acoustic
            % paths at some desired frequencies (queryFrequencies) using
            % interpolation.
            % Input arguments:
            % - acPath. (numReceiv x numSources x numFrequencies)
            % - acPathFrequencies. numFrequencies-element column vector. It
            % contains the frequencies at which the acPath has been
            % measured
            
            [numReceiv, numSources, numFreqs] = size(acPath);
            
            % Create gridded interpolant
            recPoints = 1:numReceiv;
            sourPoints = 1:numSources;
            
            gridPoints = {recPoints, sourPoints, acPathFrequencies};
            queryPoints = {recPoints, sourPoints, queryFrequencies};
            
            flag = [numReceiv, numSources, numFreqs] > 1;
            
            if numFreqs == 1
                acPath_query = repmat(acPath, [1, 1, numel(queryFrequencies)]);
            else
                gridPoints = gridPoints(flag);
                queryPoints = queryPoints(flag);
                
                F = griddedInterpolant(gridPoints, acPath, 'linear', 'nearest');
                acPath_query = F(queryPoints);
                if numReceiv == 1 && numSources > 1
                    F = griddedInterpolant(gridPoints, permute(acPath, [2 3 1]), 'linear', 'nearest');
                    acPath_query = permute(acPath_query, [3 1 2]);
                elseif numReceiv > 1 && numSources == 1
                    acPath_query = permute(acPath_query, [1 3 2]);
                end
            end
            
        end
        
        function coefficients_new = nullField(coefficients, channelMapping, frequencies, fixedChannels, acousticPathStructure)
            % Input arguments
            % - coefficients. (numSources x numFrequencies)
            
            acPath = WFSToolSimple.tuneAcousticPaths(acousticPathStructure.acousticPaths, acousticPathStructure.frequencies, frequencies); % (numReceivers x numSources x numFrequencies)
            
            % Map acoustic path channels onto coefficient channels
            [~, ind] = ismember(channelMapping, acousticPathStructure.xChannelMapping);
            acPath = acPath(:, ind, :);
            
            numFrequencies = size(coefficients, 2);
            
            fixedIndices = ismember(channelMapping, fixedChannels);
            variableIndices = ~fixedIndices; % WFS array sources
            
            coefficients_new = zeros(size(coefficients));            
            for f = 1:numFrequencies
                y = -acPath(:, fixedIndices, f)*coefficients(fixedIndices, f);
                A = acPath(:, variableIndices, f);                
                
                coefficients_new(variableIndices, f) = A\y;
            end
                      
        end
        
    end
end

