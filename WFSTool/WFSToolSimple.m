classdef WFSToolSimple < handle
    % One device, sinusoidal signals
    
    properties(SetAccess = private)
        changed % Variables changed since last update
        simplePerformance = true;
        
        % Base variables
        virtual
        real
        channelNumber
        virtualVolume
        realVolume
        signalsSpec
        sourceCorr % Correction factor from adimensional variable to physic variable
        receiverCorr % Correction factor from physics magnitude of pressure to adimiensional units
        activeReceivers % Active channels of the receiver device
        
        amplitude
        phase
        frequency
                                
        % Infrastructure variables
        fig
        player
        reprodPanel
        scenarioObj
        propPanel
        simulObj
                
        % Other
        ax
    end
    
    properties(Dependent)
        numSources
    end
    
    properties(Constant)
        c = 340; % m/s
    end
    
    % Getters and setters
    methods
        function numSources = get.numSources(obj)
            numSources = numel(obj.signalsSpec);
        end   
    end
    
    methods
        
        function obj = WFSToolSimple()
            fig = figure('Units', 'pixels', 'Position', [0 50 1200 600]);
            obj.fig = fig;
            
            obj.player = reproductorRecorder();
            obj.reprodPanel = reproductionPanel_noiseChannel(fig, @(action) obj.orderCallback(action));
            obj.scenarioObj = scenario(fig);
            obj.propPanel = propertiesPanel(fig, [0.05 0.1 0.4 0.2]);
            obj.simulObj = simulator;
            obj.ax = obj.scenarioObj.ax;
            obj.simulObj.ax = obj.ax;
            colormap(obj.ax, 'gray')
            obj.ax.CLim = [-1 1];
                        
            obj.changed = struct('virtual', false, 'real', false, 'signalsSpec', false);
            
            addlistener(obj.player, 'numChannels', 'PostSet', @(~, eventData) obj.changeScenario(eventData.AffectedObject.numChannels(1)));
            addlistener(obj.reprodPanel, 'updatedValues', @(~, evntData) obj.reprodPanelListener(evntData.type));
            addlistener(obj.player, 'playingState', 'PostSet', @(~, eventData) obj.GUIenabling(eventData.AffectedObject.playingState));
            
            obj.signalsSpec = obj.reprodPanel.signals;
            obj.virtual = obj.reprodPanel.virtual;
            obj.real = obj.reprodPanel.real;
            obj.channelNumber = obj.reprodPanel.channelNumber;
            obj.virtualVolume = obj.reprodPanel.virtualVolume;
            obj.realVolume = obj.reprodPanel.realVolume;
            obj.changed.virtual = true;
            obj.changed.real = true;
            obj.changed.signalsSpec = true;
            obj.updateEverything();
            
            obj.simulObj.XnumPoints = 200;
            obj.simulObj.YnumPoints = 200;
            obj.simulate();
            obj.simulObj.imag.HitTest = 'off';
            obj.ax.Children = flip(obj.ax.Children);
            
            uicontrol(fig, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', [0.6 0.95 0.1 0.05], 'String', 'Simulate', 'Callback', @(hObject, eventData) obj.simulateOnAxis());
        end
        
        function updateEverything(obj)
            
            if obj.changed.virtual || obj.changed.real
                rightSize = [obj.numSources, 1];
                virtualRight = all(size(obj.virtual) == rightSize);
                realRight = all(size(obj.real) == rightSize);
                signalsRight = all(size(obj.signalsSpec) == rightSize);
                channelNumberRight = all(size(obj.channelNumber) == rightSize);
                virtualVolumeRight = all(size(obj.virtualVolume) == rightSize);
                realVolumeRight = all(size(obj.virtualVolume) == rightSize);
                assert(virtualRight && realRight && signalsRight && channelNumberRight...
                    && virtualVolumeRight && realVolumeRight, 'WFSTool2:updateEverything', 'The signals specifications and the virtual and real flags must have the same size')
                
                obj.updateComMat();
            end
            
            if obj.changed.signalsSpec
                obj.updateSignalParameteres();
                obj.updateSignalProvidersVariables();
            end
        end
        
        function setNumSources(obj, numSources)
            obj.virtual = true(numSources, 1);
            obj.real = true(numSources, 1);
            obj.signalsSpec = cell(numSources, 1);
            obj.channelNumber = zeros(numSources, 1);
            obj.virtualVolume = ones(numSources, 1);
            obj.realVolume = ones(numSources, 1);
            
            for k = 1:N
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
        
        function acPath = acousticPathsTest(obj, frequency)
            
            % Generate the signal
            soundCicles = 1; % Number of sinusoidal periods for each loudspeaker
            SampleRate = 44100;
            silenceCicles = 10; % Number of periods of silence between loudspeakers
            numChannels = obj.player.numChannels(1);
            preludeSamples = SampleRate; % Prelude of samples of silence
            
            [signal, sourceCoef] = successiveChannelSinusoids( ones(numChannels, 1), -pi/2*ones(numChannels, 1), frequency, SampleRate, soundCicles, silenceCicles );
            signal = [zeros(preludeSamples, numChannels); signal];
            sourceCoef = sourceCoef*exp(-1i*2*pi*frequency*preludeSamples/SampleRate);
            
            % Reproduce and record
            obj.player.setProps('customSignal', signal, 1);
            obj.player.setProps('FsGenerator', SampleRate, 1);
            obj.player.executeOrder('play');
            
            % Analyse
            signals = obj.player.recorded(1); % Only the first recorder device (one device for reproduction, one for recording)
            SampleRateRec = obj.player.Fs_recorder(1);
            iqSignal = real2IQ(signals, SampleRateRec, frequency);
            numMicroph = size(signal, 2);
            acPath = zeros(numChannels, numMicroph);
            for l = 1:numMicroph
                % Extract the received complex coefficients
                [recCoef, ~] = pulseSignalParameters(iqSignal(:, l));
                
                % Calculate the acoustic path: relation between the source
                % coefficients and the received coefficients
                acPath(:, l) = recCoef./sourceCoef;
            end
                                   
        end
        
        function compareRecordedAndSimul(obj, frequency)
            % Real situation
            acPath = acousticPathsTest(obj, frequency);
            
            % Simulation
            microphonePos = obj.scenObj.receiversPosition;
            obj.createTheoreticalModel();
            U = calculate(obj, microphonePos); % Assume there is only one frequency (numMeasurePoints x numSources x numFrequencies)
            
            % Compare
            aux = U.'./acPath;
            
            ax = axes(figure);
            scatter(ax, 1:numChannels, abs(aux) )
        end 
        
        function fromBase2TheoreticalScenario(obj)
            
            indActiveSour = find(obj.real | obj.virtual);
            indReal = find(obj.real);
            numReal = numel(indReal);
            indRealInActive = ismember(indActiveSour, indReal);
            numLoudspeakers = obj.scenarioObj.numLoudspeakers;
           
            ldspkrsCoef = obj.getComplexCoeff();
            ldspkrsPos = obj.scenarioObj.loudspeakersPosition;
            ldspkrsOrient = simulator.vec2rotVec(obj.scenarioObj.loudspeakersOrientation);
            
            % Substitute data for the loudspeakers used as real noise
            % sources
            nSrcPos_real = obj.scenarioObj.sourcesPosition(indRealInActive, :);
            nSrcOrient_real = simulator.vec2rotVec(repmat([0 0 1], [numReal, 1]));
                    
            ldspkrsPos(obj.channelNumber(indReal), :) = nSrcPos_real;
            ldspkrsOrient(obj.channelNumber(indReal), :) = nSrcOrient_real;
            
            % Set the variables in the simulation object
            obj.simulObj.sourcePositions = ldspkrsPos;
            obj.simulObj.sourceCoefficients = ldspkrsCoef;
            obj.simulObj.sourceOrientations = ldspkrsOrient;
            obj.simulObj.radPatFuns = repmat({@(x) simulator.monopoleRadPat(x)}, [numLoudspeakers, 1]);
            obj.simulObj.freq = obj.frequency;

        end
                
        function reproduceAndRecord(obj)
            SampleRate = 44100;
            
            s = obj.getTheoreticalScenarioVariables();
            x = coefficients2signal( s.sourcesCoeff, s.frequencies, SampleRate );
                      
            obj.player.setProps('mode', 'customSignal');
            obj.player.setProps('customSignal', x, 1);
            real_comMat = obj.player.comMatrix;
            prov_comMat = false(size(real_comMat)); prov_comMat(1) = true;
            obj.player.setProps('comMatrix', prov_comMat);
            obj.player.setProps('FsGenerator', SampleRate, 1);
            obj.player.setProps('enableProc', false);
            
            obj.player.executeOrder('play'); 
            
            obj.player.setProps('enableProc', true);
            obj.player.setProps('mode', 'sinusoidal');
            obj.player.setProps('comMatrix', real_comMat);
        end
        
        function simulate(obj)
            obj.simulObj.measurePoints = obj.scenarioObj.receiversPosition;
            obj.simulObj.freq = obj.frequency;
            
            obj.simulObj.simulate();
        end
        
        function saveInformation(obj)
            defaultName = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
            [FileName,PathName, ~] = uiputfile('*.m', 'Save Information', defaultName);
            
            s = obj.exportInformation();
    
            save([PathName, FileName], s);
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
            numLoudspeakers = obj.scenarioObj.numLoudspeakers;
            numNoiseSources = obj.scenarioObj.numSources;
            numReceivers = obj.scenarioObj.numReceivers;
            
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
            s.channelReal = obj.channelNumber; % Channel for real sources
            s.freq = obj.frequency;
            s.receiverPos = obj.scenarioObj.receiversPosition;
            s.receiversOrient = simulator.vec2rotVec(repmat([0 0 1], [numReceivers, 1]));
            s.receiversRadpat = repmat({@(x) simulator.monopoleRadPat(x)}, [numReceivers, 1]); % Assumption
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
            s.receiversOrient = simulator.vec2rotVec(repmat([0 0 1], [numReceivers, 1]));
            s.receiversRadpat = repmat({@(x) simulator.monopoleRadPat(x)}, [numReceivers, 1]); % Assumption
            s.sourcesCorrectionCoeff = obj.sourceCorr;
            s.receiversCorrectionCoeff = obj.receiverCorr;   
        end
        
        function s = getExperimentalResultVariables(obj)
            recorded = obj.player.recorded{1}; % Only one recorder device
            
            s.reproducedSignal = obj.player.signalReader{1}.customSignal;
            s.reproducedSignal_SampleRate = obj.player.Fs_player(1);
            s.recordedSignal = recorded(:, obj.activeReceivers);
            s.recordedSignal_SampleRate = obj.player.Fs_recorder(1);
        end
        
        function s = getSimulationResultVariables(obj)
            s.simulatedField = obj.scenarioObj.field;
        end
    end
    
    methods(Access = private)
               
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
                    obj.virtual = obj.reprodPanel.virtual;
                    obj.real = obj.reprodPanel.real;
                    obj.signalsSpec = obj.reprodPanel.signals;
                    obj.channelNumber = obj.reprodPanel.channelNumber;
                    obj.virtualVolume = obj.reprodPanel.virtualVolume;
                    obj.realVolume = obj.reprodPanel.virtualVolume;
                    
                    obj.changed.virtual = true;
                    obj.changed.real = true;
                    obj.changed.signalsSpec = true;
                case 'virtual'
                    obj.virtual = obj.reprodPanel.virtual;
                    obj.changed.virtual = true;
                case 'real'
                    obj.real = obj.reprodPanel.real;
                    obj.changed.real = true;
                case 'channelNumber'
                    obj.channelNumber = obj.reprodPanel.channelNumber;
                    obj.updateForcedDisabledLoudspeakers();
                case 'virtualVolume'
                    obj.virtualVolume = obj.reprodPanel.virtualVolume;
                case 'realVolume'
                    obj.realVolume = obj.reprodPanel.virtualVolume;
            end
            
            obj.updateEverything();
        end
        
        function updateComMat(obj)
            % Update commutation matrix
            comMat = WFSTool_noiseChannel.createCommutationMatrix(obj.virtual, obj.real);
            obj.player.setProps('comMatrix', comMat);
            
            obj.updateSignalParameteres();
            obj.updateSignalProvidersVariables();
            obj.updateProcessorVariables_noiseChannel();
            obj.updateGUIConnectionsStuff();
            obj.updateRecorderVariables();
            
            obj.changed.real = false;
            obj.changed.virtual = false;         
        end
        
        function updateRecorderVariables(obj)
            obj.scenarioObj.setNumReceivers(numel(obj.activeReceivers));
        end
        
        function updateSignalParameteres(obj)
            
            obj.amplitude = zeros(obj.numSources, 1);
            obj.phase = zeros(obj.numSources, 1);
            obj.frequency = zeros(obj.numSources, 1);
            
            for k = 1:obj.numSources
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
            for k = 1:obj.numSources
                obj.player.setProps('mode', originType('sinusoidal'), k);
                obj.player.setProps('amplitude', obj.amplitude(k), k);
                obj.player.setProps('phase', obj.phase(k), k);
                obj.player.setProps('frequency', obj.frequency(k), k);
            end
            
            obj.changed.signalsSpec = false;
        end
                
        function updateProcessorVariables_noiseChannel(obj)
            % The processors need the getDelayFun and getAttenFun
            % functions, and these are connected to the scenario
            
            % Update scenario based on the virtual sources
            activeSources = obj.virtual | obj.real;
            obj.scenarioObj.setNumSources(sum(activeSources));
            
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
            
            obj.updateForcedDisabledLoudspeakers();
            
        end
        
        function updateForcedDisabledLoudspeakers(obj)
            numLoudspeakers = obj.player.numChannels(1);
            disabledLoudspeakers = false(numLoudspeakers, 1);
            % Assign delay and attenuation functions
            for k = 1:obj.numSources
                if obj.real(k)
                    chann = obj.channelNumber(k);
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
            
            % Update player GUI controls based on the real sources
            %             indReal = [1, find(any(obj.player.comMatrix(:, 2:end), 1))+1];
            indReal = find(any(obj.player.comMatrix, 1));
            setFrameDurationFunc = @(frameDuration) obj.player.setProps('frameDuration', frameDuration);
            setVolumeFunc = @(volume, index) obj.setVolume(volume, indReal(index));
            setDeviceFunc = @(device, index) obj.player.setProps('device', device, indReal(index));
            setDriverFunc = @(driver, index) obj.player.setProps('driver', driver, indReal(index));
            getAvailableDevicesFunc = @(index) obj.player.player{indReal(index)}.getAvailableDevices();
            labels = cell(numel(indReal), 1);
            for k = 1:numel(indReal)
                labels{k} = sprintf('Device %d', indReal(k)-1);
            end
            if ~isempty(indReal) && indReal(1) == 1
                labels{1} = 'WFS Array';
            end
            
            obj.propPanel.setFunctions(setFrameDurationFunc, setVolumeFunc, setDeviceFunc, setDriverFunc, getAvailableDevicesFunc, labels);
            
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
                virtualDelays = obj.scenarioObj.delays(:, (activeSources == activeSources(index)));
            else
                virtualDelays = zeros(numChannels, 1);
            end
            
            if obj.real(index)              
                realDelays = zeros(numChannels, 1);
                delays = [virtualDelays, realDelays];
            else
                delays = virtualDelays;
            end
        end
        
        function attenuations = getAttenuations(obj, index)
            numChannels = size(obj.scenarioObj.attenuations, 1);
            activeSources = find(obj.virtual | obj.real);
            
            if obj.virtual(index)
                virtualAttenuations = obj.virtualVolume(index)*obj.scenarioObj.attenuations(:, (activeSources == activeSources(index)));
            else
                virtualAttenuations = zeros(numChannels, 1);
            end
            
            if obj.real(index)
                realAttenuations = zeros(numChannels, 1);
                realAttenuations(obj.channelNumber(index)) = -obj.realVolume(index);
                attenuations = [virtualAttenuations, realAttenuations];
            else
                attenuations = virtualAttenuations;
            end
        end
        
        function complexCoeff = getComplexCoeff(obj)
            numFreq = obj.numSources;
            numSour = obj.scenarioObj.numLoudspeakers;
            complexCoeff = zeros(numSour, numFreq);
            for k = 1:numFreq
                delays = obj.getDelays(k);
                atten = obj.getAttenuations(k);
                
                amp = obj.amplitude(k);
                pha = obj.phase(k);
                sourceCoef = amp * exp(1i*pha);
                
                complexCoeff(:, k) = sourceCoef * atten * exp(-1i*2*pi*obj.frequency(k)*delays);
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
                        obj.fromBase2TheoreticalScenario();
                        obj.reproduceAndRecord();
                        obj.simulate();
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
        
        function changeScenario(obj, numLoudspeakers)
            numLoudspeakers = numLoudspeakers(1);
            switch numLoudspeakers
                case 0
                    sourcePosition = [0 0 0];
                    loudspeakersPosition = double.empty(0,3);
                    loudspeakersOrientation = double.empty(0,3);
                    roomPosition = [0 0 1 1];
                    obj.scenarioObj.setScenario(sourcePosition, loudspeakersPosition, loudspeakersOrientation, roomPosition);
                    
                    obj.simulObj.XLim = [roomPosition(1), roomPosition(1) + roomPosition(3)];
                    obj.simulObj.YLim = [roomPosition(2), roomPosition(2) + roomPosition(4)];
                                        
                case 2
                    sourcePosition = obj.scenarioObj.sourcesPosition;
                    loudspeakersPosition = [-0.1 0 0; 0.1 0 0];
                    loudspeakersOrientation = [1 0 0; -1 0 0];
                    roomPosition = [-2, -2, 4, 4];
                    obj.scenarioObj.setScenario(sourcePosition, loudspeakersPosition, loudspeakersOrientation, roomPosition);
                    
                    obj.simulObj.XLim = [roomPosition(1), roomPosition(1) + roomPosition(3)];
                    obj.simulObj.YLim = [roomPosition(2), roomPosition(2) + roomPosition(4)];
                    
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
                    
                    xmin = min(x); xmax = max(x); ymin = min(y); ymax = max(y);
                    xDim = xmax - xmin; yDim = ymax - ymin;
                    xmargin = 0.2 * xDim; ymargin = 0.2 * yDim;
                    roomPosition = [xmin - xmargin, ymin - ymargin, xDim + 2*xmargin, yDim + 2*ymargin];
                    
                    obj.scenarioObj.setScenario(sourcePosition, loudspeakersPosition, loudspeakersOrientation, roomPosition);
                    
                    obj.simulObj.XLim = [roomPosition(1), roomPosition(1) + roomPosition(3)];
                    obj.simulObj.YLim = [roomPosition(2), roomPosition(2) + roomPosition(4)];
                    
                otherwise
                    warning('Wrong number of output channels. There is not possible scenario for that case')
            end
            
            obj.updateForcedDisabledLoudspeakers();
        end
            
        function simulateOnAxis(obj)
            % Configure simulator object
            obj.fromBase2TheoreticalScenario();
            
            % Simulate
            obj.simulObj.generateMeasurePoints();
            obj.simulObj.simulate();
        end
        
    end
    
    methods(Static)
        function comMat = createCommutationMatrix(virtual, real)
            comMat = virtual | real;
        end
    end
end

