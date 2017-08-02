classdef WFSToolSimple < handle
    % One device, sinusoidal signals
    
    properties(SetAccess = private)
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
        WFSarrayOrientation
        WFSarrayRadiationPattern
        WFSarrayChannelMapping
        
        % Microphones
        receiverOrientation
        receiverRadiationPattern
        receiverChannelMapping % Active channels of the receiver device

        % Loudspeakers
        loudspeakerChannelMapping
        
        % Calibration
        sourceCorr % Correction factor from adimensional variable to physic variable
        receiverCorr % Correction factor from physics magnitude of pressure to adimiensional units
                       
        % Simulation results
        simulField
        
        % Experimental results
        pulseCoeffMat
        pulseLimits
        reprodFrequencies
        
        % Analysis variables
        expAcPath
        simulAcPath
        
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
        WFSarrayPosition
        
        % Microphones
        receiverPosition
        
        % Counting
        numNoiseSources
        numReceivers
        numSourcesWFSarray
        numLoudspeakers % Number of channels of the playing device
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
            numReceivers = numel(obj.receiverChannelMapping);
        end
        
        function numSourcesWFSarray = get.numSourcesWFSarray(obj)
            numSourcesWFSarray = obj.scenarioObj.numLoudspeakers;
        end
        
        function numLoudspeakers = get.numLoudspeakers(obj)
            numLoudspeakers = obj.player.numChannels(1);
        end
        
        % Noise source
        function noiseSourcePosition = get.noiseSourcePosition(obj)
            indActiveSour = find(obj.real | obj.virtual);
            indReal = find(obj.real);
            indRealInActive = ismember(indActiveSour, indReal);
            
            noiseSourcePosition = zeros(obj.numNoiseSources, 3);
            activeNoiseSourcePosition = obj.scenarioObj.sourcesPosition;
            noiseSourcePosition(indRealInActive, :) = activeNoiseSourcePosition;
        end
        
        function set.noiseSourcePosition(obj, value)
            obj.scenarioObj.sourcesPosition = value;
        end

        % WFS array
        function WFSarrayPosition = get.WFSarrayPosition(obj)
            WFSarrayPosition = obj.scenarioObj.loudspeakersPosition;
        end
        
        function set.WFSarrayPosition(obj, value)
            obj.scenarioObj.loudspeakersPosition = value;
        end

        function set.WFSarrayOrientation(obj, value)
            obj.WFSarrayOrientation = value;
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
            obj.scenarioObj.receiversPosition = value;
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
            
            
            % Noise sources
            obj.noiseSourceOrientation = [0, 0, 0, 1];
            obj.noiseSourceRadiationPattern = @(x) simulator.monopoleRadPat(x);
            
            obj.signalsSpec = obj.reprodPanel.signals;
            obj.virtual = obj.reprodPanel.virtual;
            obj.real = obj.reprodPanel.real;
            obj.noiseSourceChannelMapping = obj.reprodPanel.noiseSourceChannelMapping;
            obj.virtualVolume = obj.reprodPanel.virtualVolume;
            obj.realVolume = obj.reprodPanel.realVolume;
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
            
            for k = 1:numNoiseSources
                obj.signalsSpec{k} = '';
            end
            
            obj.updateComMat();
            obj.updateSignalProvidersVariables();
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
            sTheo = s.TheoreticalScenario;
            sExp = s.Experiment;
                                    
            % Calculate the acoustic paths
            obj.expAcPath = getAcousticPath( sTheo.frequencies, sExp.pulseCoefMat, sExp.pulseLimits, sExp.recordedSignal, sExp.recordedSignal_SampleRate);
            
            %             % Visual examination of reproduced and recorded signals
%             analyzer = WFSanalyzer();
%             analyzer.representRecordedSignal(sExp.recordedSignal, sExp.recordedSignal_SampleRate);

        end
           
        function calculateSimulatedAcousticPaths(obj)
            s = obj.exportInformation();
            sTheo = s.TheoreticalScenario;
            
            obj.simulAcPath = obj.simulObj.calculateAcousticPaths(sTheo.receiverPos);
        end
                
        function WFScalculation(obj)
                    
            ldspkrsCoef = obj.getComplexCoeff();
            ldspkrsPos = obj.WFSarrayPosition;
            ldspkrsOrient = obj.WFSarrayOrientation;
            ldspkrsRadPat = obj.WFSarrayRadiationPattern;
            
            % Substitute data for the loudspeakers used as real noise
            % sources
            indReal = find(obj.real);
            nSrcChannelMapping_real = obj.noiseSourceChannelMapping(indReal);
            nSrcCoef_real = obj.amplitude(indReal) .* exp(1i*obj.phase(indReal));
            nSrcPos_real = obj.noiseSourcePosition(indReal, :);
            nSrcOrient_real = obj.noiseSourceOrientation(indReal, :); % simulator.vec2rotVec(repmat([0 0 1], [numReal, 1]));
            nSrcRadPat_real = obj.noiseSourceRadiationPattern(indReal, :);
                             
            [flagNS, indWFS] = ismember(nSrcChannelMapping_real, obj.WFSarrayChannelMapping);
            ldspkrsCoef(indWFS, :) = nSrcCoef_real(flagNS);
            ldspkrsPos(indWFS, :) = nSrcPos_real(flagNS);
            ldspkrsOrient(indWFS, :) = nSrcOrient_real(flagNS);
            ldspkrsRadPat(indWFS) = nSrcRadPat_real(flagNS); 
            
            ldspkrsCoef = [ldspkrsCoef; nSrcCoef_real(~flagNS)];
            ldspkrsPos = [ldspkrsPos; nSrcPos_real(~flagNS)];
            ldspkrsOrient = [ldspkrsOrient; nSrcOrient_real(~flagNS)];
            ldspkrsRadPat = [ldspkrsRadPat; nSrcRadPat_real(~flagNS)];
                        
            obj.loudspeakerChannelMapping = [obj.WFSarrayChannelMapping; nSrcChannelMapping_real(~flagNS)];
            
            % Set the variables in the simulation object
            obj.simulObj.sourcePositions = ldspkrsPos;
            obj.simulObj.sourceCoefficients = ldspkrsCoef;
            obj.simulObj.sourceOrientations = ldspkrsOrient;
            obj.simulObj.radPatFuns = ldspkrsRadPat; % repmat({@(x) simulator.monopoleRadPat(x)}, [obj.numLoudspeakers, 1]);
            obj.simulObj.freq = obj.frequency;

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
            if nargin == 3
                obj.player.setProps('defaultChannelMapping', false, 1);
                obj.player.setProps('channelMapping_player', channelMapping, 1);
            end

            obj.player.executeOrder('play');
                        
            % Turn everything to the previous state
            obj.player.setProps('enableProc', true);
            obj.player.setProps('comMatrix', real_comMat);
            obj.updateSignalProvidersVariables();
            for k = 1:obj.player.numPlayers
                obj.player.setProps('driver', writingDrivers{k}, k);
                obj.player.setProps('device', writingDevices{k}, k);
            end
            
        end
        
        function reproduceAndRecord(obj)
            SampleRate = 44100;
            
            s = obj.getTheoreticalScenarioVariables();
            sourcesCoef = s.sourcesCoeff;
            frequencies = s.frequencies;
            signalFunc = @(startSample, endSample) coefficients2signal( sourcesCoef, frequencies, SampleRate, startSample, endSample);

            obj.reproduceSignalFunction(signalFunc, SampleRate, obj.loudspeakerChannelMapping);
                           
            [~, ~, obj.pulseCoeffMat, pulseLim] = coefficients2signal( sourcesCoef, frequencies, SampleRate );
            obj.pulseLimits = pulseLim/SampleRate;
        end        
        
        function simulate(obj)
            obj.simulObj.measurePoints = obj.scenarioObj.receiversPosition;
            obj.simulObj.freq = obj.frequency;
            
            obj.simulField = obj.simulObj.calculate(obj.simulObj.measurePoints);
            
        end
        
        function changeScenario(obj, numLoudspeakers)
                        
            if ~ismember(numLoudspeakers, [0 2 96])
                numLoudspeakers = 0;
            end
            
            switch numLoudspeakers
                case 0
                    sourcePosition = [0 0 0];
                    receiverPosition = zeros(obj.numReceivers, 3);
                    loudspeakersPosition = double.empty(0,3);
                    loudspeakersOrientation = double.empty(0,3);
                    roomPosition = [0 0 1 1];
                    obj.scenarioObj.setScenario(sourcePosition, receiverPosition, loudspeakersPosition, loudspeakersOrientation, roomPosition);
                    
                    obj.simulObj.XLim = [roomPosition(1), roomPosition(1) + roomPosition(3)];
                    obj.simulObj.YLim = [roomPosition(2), roomPosition(2) + roomPosition(4)];
                    delete(obj.simulObj.imag);
                case 2
                    sourcePosition = obj.scenarioObj.sourcesPosition;
                    receiverPosition = zeros(obj.numReceivers, 3);
                    loudspeakersPosition = [-0.1 0 0; 0.1 0 0];
                    loudspeakersOrientation = [1 0 0; -1 0 0];
                    roomPosition = [-2, -2, 4, 4];
                    obj.scenarioObj.setScenario(sourcePosition, receiverPosition, loudspeakersPosition, loudspeakersOrientation, roomPosition);
                    
                    obj.simulObj.XLim = [roomPosition(1), roomPosition(1) + roomPosition(3)];
                    obj.simulObj.YLim = [roomPosition(2), roomPosition(2) + roomPosition(4)];
                    delete(obj.simulObj.imag);
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
                    
                    sourcePosition = [0, 0, 0];
                    receiverPosition = zeros(obj.numReceivers, 3);
                    
                    xmin = min(x); xmax = max(x); ymin = min(y); ymax = max(y);
                    xDim = xmax - xmin; yDim = ymax - ymin;
                    xmargin = 0.2 * xDim; ymargin = 0.2 * yDim;
                    roomPosition = [xmin - xmargin, ymin - ymargin, xDim + 2*xmargin, yDim + 2*ymargin];
                    
                    obj.scenarioObj.setScenario(sourcePosition, receiverPosition, loudspeakersPosition, loudspeakersOrientation, roomPosition);
                    
                    obj.simulObj.XLim = [roomPosition(1), roomPosition(1) + roomPosition(3)];
                    obj.simulObj.YLim = [roomPosition(2), roomPosition(2) + roomPosition(4)];
                    delete(obj.simulObj.imag);
                otherwise
                    warning('Wrong number of output channels. There is not possible scenario for that case')
            end
            
            obj.updateForcedDisabledLoudspeakers();
            
            obj.sourceCorr = ones(numLoudspeakers, 1);
        end
        
        function reproduceAndRecordForAcousticPaths(obj)
            SampleRate = 44100;
            sTheo = obj.getTheoreticalScenarioVariables();
            numChannels = size(sTheo.sourcesCoeff, 1);
            frequencies = sTheo.frequencies;
            numFreq = numel(frequencies);
            
            % Set coefficients for calibration
            calCoeff = 0.1*ones(numChannels, numFreq);
            
            % Reproduce
            [ pulseCoefMat, pulseLim ] = coefficients2pulseSignalParameters( calCoeff, frequencies, SampleRate, 'preludeAndMain', 0.25, 0.25, 2 );
            signalFunc = @(startSample, endSample) pulseCoefMat2signal(frequencies, pulseCoefMat, pulseLim, SampleRate, startSample, endSample, 'sample');          
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
    end
    
    % Export and import methods
    methods
        function importInformation(obj, s)
        end
        
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
            
            % Base variables  
            sBase = obj.getBaseVariables();
            
            % Theoretical Scenario Variables
            sScen = obj.getTheoreticalScenarioVariables();

            % Simulation Results Variables
            sSimul = obj.getSimulationResultVariables();
            
            % Experimental Results Variables
            sExp = obj.getExperimentalResultVariables();
            
            s.Base = sBase;
            s.TheoreticalScenario = sScen;
            s.Simulation = sSimul;
            s.Experiment = sExp;
        end
               
        function s = getBaseVariables(obj)
            numLoudspeakers = obj.numLoudspeakers;
            numNoiseSources = obj.scenarioObj.numSources;
            
            % Base variables  
            s.loudspeakersPos = obj.scenarioObj.loudspeakersPosition;
            s.loudspeakersOrient = obj.scenarioObj.loudspeakersOrientation;
            s.loudspeakersRadPat = repmat({@(x) simulator.monopoleRadPat(x)}, [numLoudspeakers, 1]); % Assumption
            s.noiseSourcesPos = obj.scenarioObj.sourcesPosition;
            s.noiseSourcesOrient = simulator.vec2rotVec(repmat([0 0 1], [numNoiseSources, 1]));
            s.noiseSourcesRadPat = repmat({@(x) simulator.monopoleRadPat(x)}, [numLoudspeakers, 1]); % Assumption
            s.noiseSourcesCoeff = obj.amplitude.*exp(1i*obj.phase);
            s.virtual = obj.virtual; % Virtual source indices
            s.real = obj.real; % Real source indices
            s.channelReal = obj.noiseSourceChannelMapping; % Channel for real sources
            s.freq = obj.frequency;
            s.receiverPos = obj.scenarioObj.receiversPosition;
            s.receiversOrient = simulator.vec2rotVec(repmat([0 0 1], [obj.numReceivers, 1]));
            s.receiversRadpat = repmat({@(x) simulator.monopoleRadPat(x)}, [obj.numReceivers, 1]); % Assumption
            s.sourcesCorrectionCoeff = obj.sourceCorr;
            s.receiversCorrectionCoeff = obj.receiverCorr;   
            s.realVolume = obj.realVolume; 
            s.virutalVolume = obj.virtualVolume;
        end
        
        function s = getTheoreticalScenarioVariables(obj)
            s.sourcesPos = obj.simulObj.sourcePositions;
            s.sourceOrient = obj.simulObj.sourceOrientations;
            s.sourceRadPat = obj.simulObj.radPatFuns;
            s.sourcesCoeff = obj.simulObj.sourceCoefficients;
            s.frequencies = obj.frequency;
            s.receiverPos = obj.scenarioObj.receiversPosition;
            s.receiversOrient = simulator.vec2rotVec(repmat([0 0 1], [obj.numReceivers, 1]));
            s.receiversRadpat = repmat({@(x) simulator.monopoleRadPat(x)}, [obj.numReceivers, 1]); % Assumption
            s.sourcesCorrectionCoeff = obj.sourceCorr;
            s.receiversCorrectionCoeff = obj.receiverCorr;   
        end
        
        function s = getExperimentalResultVariables(obj)
            recorded = obj.player.recorded{1}; % Only one recorder device
            
            s.recordedSignal = recorded;
            s.recordedSignal_SampleRate = obj.player.Fs_recorder(1);
            s.pulseLimits = obj.pulseLimits;
            s.pulseCoefMat = obj.pulseCoeffMat;
            s.frequencies = obj.reprodFrequencies;
        end
        
        function s = getSimulationResultVariables(obj)
            s.simulatedField = obj.simulField;
        end
        
    end
    
    % WFS actions
    methods
        function WFS2realRatio(obj)
            
            realInd = obj.noiseSourceChannelMapping(1);
            WFSInd = 1:obj.simulObj.numSources;
            WFSInd(realInd) = [];
            
            U = mean(obj.simulObj.calculate(obj.simulObj.measurePoints(20000,:)), 2);
            UWFS = sum(U(WFSInd));
            Ureal = U(realInd);
            
            correction = -Ureal/UWFS;
            
            obj.simulObj.sourceCoefficients(WFSInd) = obj.simulObj.sourceCoefficients(WFSInd)*correction;
%             obj.simulObj.simulate();
            
        end
             
        function nullField(obj)
            s = obj.exportInformation();
            sTheo = s.TheoreticalScenario;
            numFreq = numel(sTheo.frequencies);
            
            sourceCoeff = sTheo.sourcesCoeff;
            
            [~, noiseOrWFSarrayFlag] = obj.mapLoudspeakersOntoNoiseAndWFSArraySources();
            indepIndices = find(noiseOrWFSarrayFlag); % Noise sources
            depIndices = find(~noiseOrWFSarrayFlag); % WFS array sources
            
            acPath = obj.expAcPath;
            sourceCoeffDep = zeros(numel(depIndices), numFreq);
            for f = 1:numFreq
                B = -acPath(indepIndices, :, f).'*sourceCoeff(indepIndices, f);
                A = acPath(depIndices, :, f).';
                
                sourceCoeffDep(:, f) = A\B;
            end
            
            newSourceCoeff = sourceCoeff;
            newSourceCoeff(depIndices) = sourceCoeffDep;
                      
            obj.simulObj.sourceCoefficients = newSourceCoeff;
        end
    end
    
    
    methods(Access = private)
               
        function timeDisplay(obj, time)
            obj.timeDisp.String = num2str(time);
        end
        
        function setDefaultValues(obj)
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
            obj.receiverChannelMapping = obj.recordPanel.activeChannels;
            obj.updateRecorderVariables();
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
                
        function updateRecorderVariables(obj)
            obj.scenarioObj.setNumReceivers(numel(obj.receiverChannelMapping));
            obj.player.setProps('channelMapping_recorder', obj.receiverChannelMapping, 1);
            obj.receiverCorr = ones(obj.numReceivers, 1);
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
            numLoudspeakers = obj.scenarioObj.numLoudspeakers;
            disabledLoudspeakers = false(numLoudspeakers, 1);
            % Assign delay and attenuation functions
            for k = 1:obj.numNoiseSources
                if obj.real(k)
                    chann = obj.noiseSourceChannelMapping(k);
                    if chann > 0 && chann <= numLoudspeakers
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
            numChannels = obj.scenarioObj.numLoudspeakers;
            activeSources = find(obj.virtual | obj.real);
            
            if obj.virtual(index)
                delays = obj.scenarioObj.delays(:, (activeSources == activeSources(index)));
            else
                delays= zeros(numChannels, 1);
            end
            
            if obj.real(index)              
                delays(obj.noiseSourceChannelMapping(index)) = 0;
            end
        end
        
        function attenuations = getAttenuations(obj, index)
            numChannels = size(obj.scenarioObj.attenuations, 1);
            activeSources = find(obj.virtual | obj.real);
            
            if obj.virtual(index)
                attenuations = obj.virtualVolume(index)*obj.scenarioObj.attenuations(:, (activeSources == activeSources(index)));
            else
                attenuations = zeros(numChannels, 1);
            end
            
            if obj.real(index)
                attenuations(obj.noiseSourceChannelMapping(index)) = -obj.realVolume(index);
            end
        end
        
        function complexCoeff = getComplexCoeff(obj)
            numFreq = obj.numNoiseSources;
            numSour = obj.scenarioObj.numLoudspeakers;
            complexCoeff = zeros(numSour, numFreq);
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
        
        function [indices, noiseOrWFSarrayFlag] = mapLoudspeakersOntoNoiseAndWFSArraySources(obj)
            % indices. numLoudspeakers-element column vector. The i-th
            % element is the corresponding theoric source (noise source or
            % WFS array/secondary source) corresponding to the i-th
            % loudspeaker
            % noiseOrWFSarrayFlag. numLoudspeakers-element column logical
            % vector. If the i-th lousdpeaker corresponds to a noise
            % source, it is set to true. If it corresponds to a noise
            % source, it is set to false.
            
            WFSarrayChanInd = (1:obj.numSourcesWFSarray)';
            noiseSourcesChanInd = obj.noiseSourceChannelMapping(obj.real);
            
            indices = zeros(obj.numLoudspeakers, 1);
            indices(WFSarrayChanInd) = (1:numel(WFSarrayChanInd))';
            indices(noiseSourcesChanInd) = (1:numel(noiseSourcesChanInd))';
            
            noiseOrWFSarrayFlag = false(size(indices));
            noiseOrWFSarrayFlag(noiseSourcesChanInd) = true;
            
        end
           
        function simulateOnAxis(obj)
            % Configure simulator object
            obj.WFScalculation();
            
            % Simulate
            obj.simulObj.generateMeasurePoints();
            obj.simulObj.simulate();
        end
        
    end
    
    methods(Static)
        function comMat = createCommutationMatrix(virtual, ~)
%             comMat = virtual | real;
            comMat = true(numel(virtual), 1);
        end
    end
end

