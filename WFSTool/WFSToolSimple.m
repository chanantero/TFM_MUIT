classdef WFSToolSimple < handle
    % One device, sinusoidal signals
    
    properties(SetAccess = public)
        changed % Variables changed since last update
        simplePerformance = true;
        
        % Noise source
        virtual
        virtualVolume
        real
        realVolume
        noiseSourceChannelMapping
        signalsSpec % String that specifies the coefficient
        amplitude
        phase
        
        % WFS Array
        WFSarrayChannelMapping
        
        % Microphones
        receiverChannelMapping % Active channels of the receiver device

        % Loudspeakers
        loudspeakerChannelMapping
        loudspeakerMapping % Structure
                                   
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
        noiseSourceAcPathStruct
        WFSarrayAcPathStruct

        % Infrastructure variables
        fig
        player
        reprodPanel
        scenarioObj
        propPanel
        recordPanel
        simulTheo
                
        % Other
        ax
        timeDisp
    end
    
    properties(Dependent)
        % Noise source
        noiseSourcePosition
        noiseSourceCoefficient % Without specifying the frequencies. (obj.numNoiseSources x 1)
        noiseSourceCoefficient_complete % All frequencies. (obj.numNoiseSources x obj.numNoiseSources)
        noiseSourceOrientation
        noiseSourceRadiationPattern
        frequency
        noiseSourceIndSimulTheo
   
        % WFS Array
        WFSarrayPosition
        WFSarrayOrientation
        WFSarrayRadiationPattern
        WFSarrayCoefficient
        WFSarrayIndSimulTheo
        
        % Microphones
        receiverPosition
        receiverOrientation
        receiverRadiationPattern
        
        % Loudspeakers
        loudspeakerCoefficient
        
        % Experiment
        recorded
        
        % Simulation
        simulField
        
        % Counting
        numNoiseSources
        numReceivers
        numSourcesWFSarray
        numSourcesTheo
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
            numNoiseSources = numel(obj.noiseSourceChannelMapping);
        end   
        
        function numReceivers = get.numReceivers(obj)
            numReceivers = numel(obj.receiverChannelMapping);
        end
        
        function numSourcesWFSarray = get.numSourcesWFSarray(obj)
            numSourcesWFSarray = numel(obj.WFSarrayChannelMapping);
        end
        
        function numLoudspeakers = get.numLoudspeakers(obj)
            numLoudspeakers = numel(obj.loudspeakerChannelMapping);
        end
        
        function numSourcesTheo = get.numSourcesTheo(obj)
            numSourcesTheo = obj.numNoiseSources + obj.numSourcesWFSarray;
        end
        
        
        % Noise source
        function noiseSourcePosition = get.noiseSourcePosition(obj)
%             indActiveSour = find(obj.real | obj.virtual);
% %             indReal = find(obj.real);
% %             indRealInActive = ismember(indActiveSour, indReal);
%             
%             noiseSourcePosition = zeros(obj.numNoiseSources, 3);
%             noiseSourcePosition(indActiveSour, :) = obj.scenarioObj.sourcesPosition(indActiveSour, :);
            noiseSourcePosition = obj.simulTheo.sourcePositions(obj.noiseSourceIndSimulTheo, :);
        end
        
        function set.noiseSourcePosition(obj, value)
            numPos = size(value, 1);
            
            if numPos == obj.numNoiseSources
                obj.scenarioObj.setSourcePosition(value, (1:size(value, 1)));
            else
                warning('WFSToolSimple:WrongNumberNoiseSources', 'The number of positions must be equal to the number of sources');
            end
            
            % Set it in the simulator object
            obj.simulTheo.sourcePositions(obj.noiseSourceIndSimulTheo, :) = value;
        end
        
        function noiseSourceCoefficient = get.noiseSourceCoefficient(obj)
            noiseSourceCoefficient = obj.amplitude .* exp(1i*obj.phase);
        end
        
        function set.noiseSourceCoefficient(obj, value)
            obj.amplitude = abs(value);
            obj.phase = angle(value);
        end
        
        function noiseSourceCoefficient_complete = get.noiseSourceCoefficient_complete(obj)
            noiseSourceCoefficient_complete = obj.simulTheo.sourceCoefficients(obj.noiseSourceIndSimulTheo, :);
        end
        
        function set.noiseSourceCoefficient_complete(obj, value)
            obj.simulTheo.sourceCoefficients(obj.noiseSourceIndSimulTheo, :) = value;
        end
        
        function noiseSourceOrientation = get.noiseSourceOrientation(obj)
            noiseSourceOrientation = obj.simulTheo.sourceOrientations(obj.noiseSourceIndSimulTheo, :);
        end
        
        function set.noiseSourceOrientation(obj, value)
            obj.simulTheo.sourceOrientations(obj.noiseSourceIndSimulTheo, :) = value;
        end
        
        function noiseSourceRadiationPattern = get.noiseSourceRadiationPattern(obj)
            noiseSourceRadiationPattern = obj.simulTheo.radPatFuns(obj.noiseSourceIndSimulTheo);
        end
        
        function set.noiseSourceRadiationPattern(obj, value)
            obj.simulTheo.radPatFuns(obj.noiseSourceIndSimulTheo) = value;
        end
        
        function frequency = get.frequency(obj)
            frequency = obj.simulTheo.freq;
        end
        
        function set.frequency(obj, value)
            obj.simulTheo.freq = value;
        end

        function noiseSourceIndSimulTheo = get.noiseSourceIndSimulTheo(obj)
            noiseSourceIndSimulTheo = obj.numSourcesWFSarray + (1:obj.numNoiseSources);
        end
        

        % WFS array
        function WFSarrayPosition = get.WFSarrayPosition(obj)
            WFSarrayPosition = obj.simulTheo.sourcePositions(obj.WFSarrayIndSimulTheo, :);
        end
        
        function set.WFSarrayPosition(obj, value)
            obj.scenarioObj.setLoudspeakerPosition(value, 1:obj.numSourcesWFSarray);
            obj.simulTheo.sourcePositions(obj.WFSarrayIndSimulTheo, :) = value;
        end

        function set.WFSarrayOrientation(obj, value)
            obj.scenarioObj.loudspeakersOrientation = simulator.rotVec2BroadsideVec(value);
            obj.simulTheo.sourceOrientations(obj.WFSarrayIndSimulTheo, :) = value;
        end

        function WFSarrayOrientation = get.WFSarrayOrientation(obj)
            WFSarrayOrientation = obj.simulTheo.sourceOrientations(obj.WFSarrayIndSimulTheo, :);
%             WFSarrayOrientation = simulator.vec2rotVec(obj.scenarioObj.loudspeakersOrientation);
        end
        
        function WFSarrayRadiationPattern = get.WFSarrayRadiationPattern(obj)
            WFSarrayRadiationPattern = obj.simulTheo.radPatFuns(obj.WFSarrayIndSimulTheo);
        end
        
        function set.WFSarrayRadiationPattern(obj, value)
            obj.simulTheo.radPatFuns(obj.WFSarrayIndSimulTheo) = value;
        end
        
        function WFSarrayCoefficient = get.WFSarrayCoefficient(obj)
            WFSarrayCoefficient = obj.simulTheo.sourceCoefficients(obj.WFSarrayIndSimulTheo, :);
        end
        
        function set.WFSarrayCoefficient(obj, value)          
            obj.simulTheo.sourceCoefficients(obj.WFSarrayIndSimulTheo, :) = value;
        end
        
        function WFSarrayIndSimulTheo = get.WFSarrayIndSimulTheo(obj)
            WFSarrayIndSimulTheo = 1:obj.numSourcesWFSarray;
        end
        
        % Microphones
        function receiverPosition = get.receiverPosition(obj)
            receiverPosition = obj.scenarioObj.receiversPosition;
        end
        
        function set.receiverPosition(obj, value)
            obj.scenarioObj.setReceiverPosition(value, (1:size(value, 1)));
            obj.simulTheo.recPositions = value;
        end
           
        function receiverOrientation = get.receiverOrientation(obj)
            receiverOrientation = obj.simulTheo.recOrientations;
        end
        
        function set.receiverOrientation(obj, value)
            obj.simulTheo.recOrientations = value;
        end
         
        function receiverRadiationPattern = get.receiverRadiationPattern(obj)
            receiverRadiationPattern = obj.simulTheo.recRadPatFuns;
        end
        
        function set.receiverRadiationPattern(obj, value)
            obj.simulTheo.recRadPatFuns = value;
        end
        
         
        % Loudspeakers
        function loudspeakerCoefficient = get.loudspeakerCoefficient(obj)
            % Calculate the loudspeaker coefficients
            mapping = obj.loudspeakerMapping;
            
            loudspeakerCoefficient = zeros(obj.numLoudspeakers, obj.numNoiseSources);
            loudspeakerCoefficient(mapping(1).destinationInd, :) = obj.WFSarrayCoefficient(mapping(1).originInd, :);
            loudspeakerCoefficient(mapping(2).destinationInd, :) = obj.noiseSourceCoefficient_complete(mapping(2).originInd, :);
        end
                      
        % Experimental
        function recorded = get.recorded(obj)
            recorded = obj.player.recorded{1};
        end
                
        % Simulation
        function simulField = get.simulField(obj)
            simulField = obj.simulTheo.field;
        end
        
        % Channel mappings
        function set.noiseSourceChannelMapping(obj, value)
            obj.noiseSourceChannelMapping = value;
            obj.updateLoudspeakerMappingVariables();
        end
        
        function set.WFSarrayChannelMapping(obj, value)
            obj.WFSarrayChannelMapping = value;
            obj.updateLoudspeakerMappingVariables();
        end
    end
    
    % Export methods
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
            s.loudspeakerCoefficient = obj.loudspeakerCoefficient;
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
            obj.simulTheo = simulator;
            obj.ax = obj.scenarioObj.ax;
            obj.simulTheo.ax = obj.ax;
            colormap(obj.ax, 'gray')
            obj.ax.CLim = [-1 1];            
                        
            obj.changed = struct('virtual', false, 'real', false, 'signalsSpec', false);
            
            addlistener(obj.reprodPanel, 'updatedValues', @(~, evntData) obj.reprodPanelListener(evntData.type));
            addlistener(obj.recordPanel, 'updatedValues', @(~, evntData) obj.recordPanelListener(evntData.type));
            addlistener(obj.scenarioObj, 'updatedValues', @(~, ~) obj.updateVariablesBasedOnScenario());
            addlistener(obj.player, 'playingState', 'PostSet', @(~, eventData) obj.GUIenabling(eventData.AffectedObject.playingState));
            addlistener(obj.player, 'numChannels', 'PostSet', @(~, eventData) obj.setNumWFSarraySources(eventData.AffectedObject.numChannels(1)));
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
            
            obj.simulTheo.XnumPoints = 200;
            obj.simulTheo.YnumPoints = 200;
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
            obj.adjustSimulTheoSources(obj.numSourcesWFSarray, obj.numNoiseSources, 'noiseSource', numNoiseSources);
            obj.noiseSourceChannelMapping = zeros(numNoiseSources, 1);
            obj.updateLoudspeakerMappingVariables();
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
        
        function setNumWFSarraySources(obj, numWFSsources)
                  
            s = WFSToolSimple.generateScenario(numWFSsources);
            obj.adjustSimulTheoSources(obj.numSourcesWFSarray, obj.numNoiseSources, 'WFSarray', numWFSsources);
            obj.WFSarrayChannelMapping = (1:numWFSsources)';
            obj.updateLoudspeakerMappingVariables();
            obj.scenarioObj.setScenario(obj.noiseSourcePosition, obj.receiverPosition, s.loudspeakersPosition, s.loudspeakersOrientation, s.roomPosition);
            obj.WFSarrayPosition = s.loudspeakersPosition;
            obj.WFSarrayRadiationPattern = repmat({@ simulator.monopoleRadPat}, numWFSsources, 1);
            uniqueFreq = unique(obj.frequency);
            obj.WFSarrayAcPathStruct = struct('acousticPaths', zeros(obj.numReceivers, numWFSsources, numel(uniqueFreq)), 'frequencies', uniqueFreq);
            
            
            obj.simulTheo.XLim = [s.roomPosition(1), s.roomPosition(1) + s.roomPosition(3)];
            obj.simulTheo.YLim = [s.roomPosition(2), s.roomPosition(2) + s.roomPosition(4)];
            delete(obj.simulTheo.imag);
            
            obj.updateForcedDisabledLoudspeakers();
            
            obj.WFSarrayCoefficient = zeros(numWFSsources, obj.numNoiseSources);
            
        end
        
        function setNumReceivers(obj, numReceivers)
            obj.scenarioObj.setNumReceivers(numReceivers);
            obj.receiverPosition = obj.receiverPosition;
            obj.receiverOrientation = repmat([1, 0, 0, 1], numReceivers, 1);
            obj.receiverRadiationPattern = repmat({@ simulator.monopoleRadPat}, numReceivers);
            obj.receiverChannelMapping = (1:numReceivers)';
            obj.player.setProps('channelMapping_recorder', obj.receiverChannelMapping, 1);
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
            obj.updateLoudspeakerMappingVariables();
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
            
            % Calculate coefficients
                % A) Set noise source coefficients
                obj.noiseSourceCoefficient_complete = diag(obj.noiseSourceCoefficient);

                % B) Set WFS array coefficients
                obj.WFSarrayCoefficient = obj.getComplexCoeffWFS();

                % C) Adjust WFS array coefficients
                obj.tuneWFSField();

                % D) Apply real and virtual flags
                obj.noiseSourceCoefficient_complete(:, ~obj.real) = 0;
                obj.WFSarrayCoefficient(:, ~obj.virtual) = 0;
            
            % Set acoustic paths
            WFSarrayAcPath = WFSToolSimple.tuneAcousticPaths(obj.WFSarrayAcPathStruct.acousticPaths, obj.WFSarrayAcPathStruct.frequencies, obj.frequency);
            noiseSourcesAcPath = WFSToolSimple.tuneAcousticPaths(obj.noiseSourceAcPathStruct.acousticPaths, obj.noiseSourceAcPathStruct.frequencies, obj.frequency);
            
            acPath = zeros(obj.numReceivers, obj.numSourcesTheo, obj.numNoiseSources);
            acPath(:, obj.WFSarrayIndSimulTheo, :) = WFSarrayAcPath;
            acPath(:, obj.noiseSourceIndSimulTheo, :) = noiseSourcesAcPath;
            
            obj.simulTheo.acPath = acPath;
                     
        end
        
        function updateLoudspeakerMappingVariables(obj)
            indReal = find(obj.real);
            noiseSourceChannelMapping_real = obj.noiseSourceChannelMapping(obj.real);
            [obj.loudspeakerChannelMapping, mapping] = combineIndices(obj.WFSarrayChannelMapping, noiseSourceChannelMapping_real);
            mapping(2).originInd = indReal(mapping(2).originInd);
            obj.loudspeakerMapping = mapping;
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
            obj.player.setProps('channelMapping_recorder', obj.receiverChannelMapping, 1);
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
            frequencies = obj.frequency;
            
            [ pulseCoefMat, pulseLim ] = coefficients2pulseSignalParameters( obj.loudspeakerCoefficient, frequencies, SampleRate, option, varargin{:} );
            signalFunc = @(startSample, endSample) pulseCoefMat2signal(pulseCoefMat, pulseLim, frequencies, SampleRate, startSample, endSample);           
%             signalFunc = @(startSample, endSample) coefficients2signal( sourcesCoef, frequencies, SampleRate, startSample, endSample, false, option, varargin{:});

            obj.reproduceSignalFunction(signalFunc, SampleRate, obj.loudspeakerChannelMapping);
                           
            obj.pulseCoeffMat = pulseCoefMat;
            obj.pulseLimits = pulseLim/SampleRate;
            obj.reprodFrequencies = frequencies;
            
        end        
        
        function simulate(obj)
            % Simulate
            obj.simulTheo.updateField();
            
%             % Draw
%             obj.simulTheo.drawField();           
        end        
        
        function reproduceAndRecordForAcousticPaths(obj)
            SampleRate = 44100;
            numChannels = obj.numLoudspeakers;
            frequencies = obj.frequency;
            uniqueFreq = unique(frequencies);
            numUniqueFreq = numel(uniqueFreq);
            
            % Set coefficients for calibration          
            calCoeff = 0.4*ones(numChannels, numUniqueFreq);
            
            % Reproduce
            [ pulseCoefMat, pulseLim ] = coefficients2pulseSignalParameters( calCoeff, uniqueFreq, SampleRate, 'preludeAndMain', 1, 1, 2 );
            signalFunc = @(startSample, endSample) pulseCoefMat2signal(pulseCoefMat, pulseLim, uniqueFreq, SampleRate, startSample, endSample, 'type_marker', 'sample');          
            obj.reproduceSignalFunction(signalFunc, SampleRate);
            
            % Save information about the reproduced signal
            obj.pulseCoeffMat = pulseCoefMat;
            obj.pulseLimits = pulseLim/SampleRate;
            obj.reprodFrequencies = uniqueFreq;
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
            obj.updateForcedDisabledLoudspeakers();
        end
        
        function updateRecordPanelBasedOnVariables(obj)
           obj.recordPanel.setActiveChannels(obj.receiverChannelMapping);
        end
        
        function theoricWFSacousticPath(obj)
            
            uniqueFreq = unique(obj.frequency);
            
            acPath = simulator.calculateTheoricAcousticPaths(...
                obj.WFSarrayPosition, obj.WFSarrayRadiationPattern, obj.WFSarrayOrientation,...
                obj.receiverPosition, obj.receiverRadiationPattern, obj.receiverOrientation, uniqueFreq, obj.c);
            
            obj.WFSarrayAcPathStruct = struct('acousticPaths', acPath, 'frequencies', uniqueFreq);
            
        end
        
        function theoricNoiseSourceAcousticPath(obj)
            
            uniqueFreq = unique(obj.frequency);

            acPath = simulator.calculateTheoricAcousticPaths(...
                obj.noiseSourcePosition, obj.noiseSourceRadiationPattern, obj.noiseSourceOrientation,...
                obj.receiverPosition, obj.receiverRadiationPattern, obj.receiverOrientation, uniqueFreq, obj.c);
            
            obj.noiseSourceAcPathStruct = struct('acousticPaths', acPath, 'frequencies', uniqueFreq);
            
        end
        
    end
    
    % WFS actions
    methods
             
        function tuneWFSField(obj)
            % When the coefficients returned by the scenario object don't cancel the field,
            % this function can multiply all of them by a single complex value in order to create a
            % field that cancels the original field.
            
            % This is only applied to components of the loudspeaker
            % coefficients that come from virtual noise sources, as the WFS
            % processing only is relevant on them.
            
            % Define a set of points where the theoric field should be
            % null because of the cancellation
            xLim = [1, 2]; yLim = [2, 4];
            x = linspace(xLim(1), xLim(2), 10);
            y = linspace(yLim(1), yLim(2), 10);
            z = 0;
            [X, Y, Z] = ndgrid(x, y, z);
            testPoints = [X(:), Y(:), Z(:)];
            
            % Calculate the contribution of each source to each receiver
            % point
            % Use a simulator object of the theoric scenario. In this
            % object, all noise sources are situated, as well as all WFS
            % array sources. All noise sources are active.
%             U = calculateField(acousticPaths, sourceCoefficients, true);
            
            U = obj.simulTheo.calculateTheoricField(testPoints, true); % (numTestPoints x obj.numSourcesTheo x obj.numNoiseSources)            
            
            % Apply the correction/tunning for each layer on the 3rd
            % dimension, i.e., for each noise source that is virtual
            
            % Find which of the loudspeakers correspond to real sources and
            % which of them belong to the WFS array
            realFlag = obj.noiseSourceIndSimulTheo;
            WFSflag = obj.WFSarrayIndSimulTheo;
            
            numVirtualSources = sum(obj.virtual);
            U_WFS = sum(U(:, WFSflag, obj.virtual), 2); % Contribution of all WFS array sources to each test point. (numTestPoints x 1 x obj.numNoiseSources)
            U_real = sum(U(:, realFlag, obj.virtual), 2); % Contribution of all noise sources. (numTestPoints x 1 x obj.numNoiseSources)
            
            % Find a coefficient C that minimizes the squared error: (U_WFS * C + U_real).^2
            C = zeros(1, numVirtualSources, 1);
            for k = 1:numVirtualSources
                C(k) = -U_WFS(:, 1, k)\U_real(:, 1, k);
            end
                  
            obj.WFSarrayCoefficient(:, obj.virtual) = obj.WFSarrayCoefficient(:, obj.virtual).*repmat(C, [obj.numSourcesWFSarray, 1]);
        end
        
        function cancellField(obj)
                                  
            coefficients_new = WFSToolSimple.nullField(obj.loudspeakerCoefficient,...
                obj.loudspeakerChannelMapping, obj.loudspeakerFrequencies,...
                obj.noiseSourceChannelMapping(obj.real), obj.loudspeakerAcPathStruct);
            
            obj.loudspeakerCoefficient = coefficients_new;
            
        end
                   
    end    
    
    methods(Access = private)
        
        function adjustSimulTheoSources(obj, numWFSSour, numNoiseSour, newType, numSour)
            switch newType
                case 'noiseSource'
                    indWFS = 1:numWFSSour;
                    WFSPos = obj.simulTheo.sourcePositions(indWFS, :);
                    WFSOrient = obj.simulTheo.sourceOrientations(indWFS, :);
                    WFSRadPat = obj.simulTheo.radPatFuns(indWFS, :);
                    
                    obj.simulTheo.sourcePositions = [WFSPos; zeros(numSour, 3)];
                    obj.simulTheo.sourceCoefficients = zeros(numSour + numWFSSour, numSour); % Reset
                    obj.simulTheo.sourceOrientations = [WFSOrient; zeros(numSour, 4)];
                    obj.simulTheo.radPatFuns = [WFSRadPat; cell(numSour, 1)];

                case 'WFSarray'
                    indNoise = numWFSSour + (1:numNoiseSour);
                    noisePos = obj.simulTheo.sourcePositions(indNoise, :);
                    noiseCoef = obj.simulTheo.sourceCoefficients(indNoise, :);
                    noiseOrient = obj.simulTheo.sourceOrientations(indNoise, :);
                    noiseRadPat = obj.simulTheo.radPatFuns(indNoise, :);
                    
                    obj.simulTheo.sourcePositions = [zeros(numSour, 3); noisePos];
                    obj.simulTheo.sourceCoefficients = [zeros(numSour, numNoiseSour); noiseCoef];
                    obj.simulTheo.sourceOrientations = [zeros(numSour, 4); noiseOrient];
                    obj.simulTheo.radPatFuns = [cell(numSour, 1); noiseRadPat];
            end
            
        end
        
        function updateVariablesBasedOnScenario(obj)
            activeSources = obj.virtual | obj.real; % This are the sources that are in the scenario object
            obj.noiseSourcePosition(activeSources, :) = obj.scenarioObj.sourcesPosition;
            obj.receiverPosition = obj.scenarioObj.receiversPosition;
        end
        
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
            freq = zeros(obj.numNoiseSources, 1);
            
            for k = 1:obj.numNoiseSources
                signalSpec = obj.signalsSpec{k};
                
                % Is a complex number?
                param = regexp(signalSpec, 'A:(?<Amplitude>(\d+\.\d+|\d+)) Ph:(?<Phase>(\d+\.\d+|\d+)) f:(?<Frequency>(\d+\.\d+|\d+))', 'names');
                if isempty(param)
                    % It is not a complex number
                    % As default, treat it as a tone
                    obj.amplitude(k) = 0;
                    obj.phase(k) = 0;
                    freq(k) = 1;
                else
                    % It is a complex number, set the
                    % properties
                    obj.amplitude(k) = str2double(param.Amplitude);
                    obj.phase(k) = str2double(param.Phase);
                    freq(k) = str2double(param.Frequency);
                end
            end
            
            obj.frequency = freq;
            
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
%             delayShift = angle(obj.sourceCorr)/(2*pi*obj.frequency(index));
%             delays = delays + delayShift;
            
        end
        
        function attenuations = attenuationFunction(obj, index)
            
            attenuations = obj.getAttenuations(index); %./abs(obj.sourceCorr);
            
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
        
        function delays = getDelaysWFS(obj, indices)
                        
            activeSources = find(obj.virtual | obj.real);
            [isActive, indScenario] = ismember(indices, activeSources);
            
            delays = zeros(obj.numSourcesWFSarray, numel(indices));
            delays(:, isActive) = obj.scenarioObj.delays(:, indScenario(isActive));
                        
%             if obj.virtual(index)
%                 % Find index of source in scenarioObj. Only the sources that
%                 % are active (virtual, real or both flags set to true) are in
%                 % scenarioObj
%                 activeSources = find(obj.virtual | obj.real);
%                 ind = (activeSources == activeSources(index));
%                 
%                 delays = obj.scenarioObj.delays(:, ind);
%             else
%                 delays= zeros(obj.numSourcesWFSarray, 1);
%             end
            
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
        
        function attenuations = getAttenuationsWFS(obj, indices)
                    
            activeSources = find(obj.virtual | obj.real);
            [isActive, indScenario] = ismember(indices, activeSources);
            
            attenuations = zeros(obj.numSourcesWFSarray, numel(indices));
            attenuations(:, isActive) = obj.scenarioObj.attenuations(:, indScenario(isActive));
                                    
        end
        
        function complexCoeff = getComplexCoeffWFS(obj)
            complexCoeff = zeros(obj.numSourcesWFSarray, obj.numNoiseSources);
            
            for k = 1:obj.numNoiseSources
                delays = obj.getDelaysWFS(k);
                atten = obj.getAttenuationsWFS(k);                                
                sourceCoef = obj.noiseSourceCoefficient(k);
                
                complexCoeff(:, k) = sourceCoef * atten .* exp(-1i*2*pi*obj.frequency(k)*delays);
            end          
        end
               
        function complexCoeff = getComplexCoeff(obj)
            complexCoeff = zeros(obj.numLoudspeakers, obj.numNoiseSources);
            for k = 1:obj.numNoiseSources
                delays = obj.getDelays(k);
                atten = obj.getAttenuations(k);
                sourceCoef = obj.noiseSourceCoefficient(k);
                
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
            obj.simulTheo.generateMeasurePoints();
            obj.simulTheo.updateTheoricAcousticPathsImage();
            obj.simulTheo.updateFieldImage();
            
            % Draw
            obj.simulTheo.drawImage();
        end
        
    end
    
    % Static methods
    methods(Static)
        
        function comMat = createCommutationMatrix(virtual, ~)
%             comMat = virtual | real;
            comMat = true(numel(virtual), 1);
        end
        
        function s = generateScenario(numLoudspeakers)
                        
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
                    loudspeakersPosition = zeros(numLoudspeakers, 3);
                    loudspeakersOrientation = repmat([1 0 0], [numLoudspeakers, 1]);
                    roomPosition = [-1, -1, 2, 2];
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
                    acPath_query = permute(acPath_query, [3 1 2]);
                elseif numReceiv > 1 && numSources == 1
                    acPath_query = permute(acPath_query, [1 3 2]);
                end
            end
            
        end
        
        function coefficients_new = nullField(coefficients, channelMapping, frequencies, fixedChannels, acousticPathStructure)
            % Input arguments
            % - coefficients. (numSources x numFrequencies)
            % - channelMapping. (numSources x 1)
            % - frequencies. (numFrequencies x 1)
            % - fixedChannels. (numFixedChannels x 1)
            % - acousticPathStructure.
            
            numSources = size(coefficients, 1);
            numFrequencies = size(coefficients, 2);
            
            acPath = WFSToolSimple.tuneAcousticPaths(acousticPathStructure.acousticPaths, acousticPathStructure.frequencies, frequencies); % (numReceivers x numSources x numFrequencies)
            numRec = size(acPath, 1);
            
            % Map acoustic path channels onto coefficient channels
            if isfield(acousticPathStructure, 'xChannelMapping')
                acPathXchannelMapping = acousticPathStructure.xChannelMapping;
            else
                acPathXchannelMapping = 1:size(acPath, 2);
            end
            [flag, ind] = ismember(channelMapping, acPathXchannelMapping);
            
            acPath_adapt = zeros(numRec, numSources, numFrequencies);
            acPath_adapt(:, flag, :) = acPath(:, ind(flag), :);
            
            % Set fixed and variable flags
            fixedIndices = ismember(channelMapping, fixedChannels);
            variableIndices = ~fixedIndices; % WFS array sources
            
            coefficients_new = zeros(size(coefficients));            
            for f = 1:numFrequencies
                y = -acPath_adapt(:, fixedIndices, f)*coefficients(fixedIndices, f);
                A = acPath_adapt(:, variableIndices, f);                
                
                coefficients_new(variableIndices, f) = A\y;
                coefficients_new(fixedIndices, f) = coefficients(fixedIndices, f);
            end
                      
        end
        
    end
end

