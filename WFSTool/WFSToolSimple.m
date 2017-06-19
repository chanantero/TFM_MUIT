classdef WFSToolSimple < handle
    % One device, sinusoidal signals
    
    properties(SetAccess = private)
        changed % Variables changed since last update
        
        % High level variables
        virtual
        real
        channelNumber
        virtualVolume
        realVolume
                
        % Scenario variables
        signalsSpec
        amplitude
        phase
        frequency
        sourceCorr % Correction factor from adimensional variable to physic variable
        receiverCorr % Correction factor from physics magnitude of pressure to adimiensional units
                                
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
        numReceivers
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
            
            obj.simulObj.XnumPoints = 100;
            obj.simulObj.YnumPoints = 100;
            obj.simulate();
            obj.simulObj.imag.HitTest = 'off';
            obj.ax.Children = flip(obj.ax.Children);
            
            uicontrol(fig, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', [0.6 0.95 0.1 0.05], 'String', 'Simulate', 'Callback', @(hObject, eventData) obj.simulate());
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
        
        function acPath = acousticPathsTest(obj, frequency)
            
            % Generate the signal
            soundCicles = 1; % Number of sinusoidal periods for each loudspeaker
            SampleRate = 44100;
            silenceCicles = 10; % Number of periods of silence between loudspeakers
            numChannels = obj.player.numChannels(1);
            preludeSamples = SampleRate; % Prelude of samples of silence
            
            cicleSamples = ceil(SampleRate/frequency*soundCicles);
            silenceSamples = ceil(SampleRate/frequency*silenceCicles);
            samplesPerChannel = cicleSamples + silenceSamples;
            
            signal = zeros(preludeSamples + numChannels*samplesPerChannel, numChannels);
            for k = 1:numChannels
                ind = preludeSamples + (k - 1)*samplesPerChannel + (1:cicleSamples);
                t = (0:cicleSamples - 1)/SampleRate;
                signal(ind, k) = sin(2*pi*frequency*t);
            end
            
            % Determine the complex coefficients of the source signal for
            % each channel
            amp = onex(numChannels, 1);
            t_ini = (preludeSamples + samplesPerChannel*(0:numChannels - 1))/SampleRate;
            phas = 2*pi*frequency*t_ini - pi/2;
            sourceCoef = amp*exp(1i*phas);
            
            % Reproduce and record
            obj.player.setProps('customSignal', signal, 1);
            obj.player.setProps('FsGenerator', SampleRate, 1);
            obj.player.executeOrder('play');
            
            % Analyse
            signals = obj.player.recorded;
            SampleRateRec = obj.player.Fs_recorder;
            numRecorders = numel(signals);
            for k = 1:numRecorders
                iqSignal = real2IQ(signals{k}, SampleRateRec(k), frequency);
                
                % Extract the received complex coefficients
                [recCoef, ~] = pulseSignalParameters(iqSignal);
                
                % Calculate the acoustic path: relation between the source
                % coefficients and the received coefficients
                acPath = recCoef./sourceCoef;
            end
                        
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
            
            obj.changed.real = false;
            obj.changed.virtual = false;         
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
            numChannels = obj.scenarioObj.numLoudspeakers;
            activeSources = find(obj.virtual | obj.real);
            
            if obj.virtual(index)
                virtualDelays = obj.scenarioObj.delays(:, (activeSources == activeSources(index)));
            else
                virtualDelays = zeros(numChannels, 1);
            end
            
            if obj.real(index)
                % delays(obj.channelNumber(index)) = 0;
                realDelays = zeros(numChannels, 1);
                delays = [virtualDelays, realDelays];
            else
                delays = virtualDelays;
            end
            
        end
        
        function attenuations = attenuationFunction(obj, index)
            numChannels = size(obj.scenarioObj.attenuations, 1);
            activeSources = find(obj.virtual | obj.real);
            
            if obj.virtual(index)
                virtualAttenuations = obj.virtualVolume(index)*obj.scenarioObj.attenuations(:, (activeSources == activeSources(index)));
            else
                virtualAttenuations = zeros(numChannels, 1);
            end
            
            if obj.real(index)
                % attenuations(obj.channelNumber(index)) = obj.realVolume(index);
                realAttenuations = zeros(numChannels, 1);
                realAttenuations(obj.channelNumber(index)) = -obj.realVolume(index);
                attenuations = [virtualAttenuations, realAttenuations];
            else
                attenuations = virtualAttenuations;
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
                    obj.player.executeOrder('play');
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
        
        function changeScenario(obj, numChannels)
            numChannels = numChannels(1);
            switch numChannels
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
              
        function createTheoreticalModel(obj)
             
            indActiveSour = find(obj.real | obj.virtual);
            numActiveSour = numel(indActiveSour);
            indReal = find(obj.real);
            numReal = numel(indReal);
            indRealInActive = ismember(indActiveSour, indReal);
            indVirtual = find(obj.virtual);
            numVirtual = numel(indVirtual);
            indVirtualInActive = ismember(indActiveSour, indVirtual);
            
            % Active sources
            ampReal = obj.amplitude(indReal).*obj.realVolume(indReal);
            phaseReal = obj.phase(indReal);
            nSrcCoef_real = ampReal.*exp(1i*phaseReal);
            
            nSrcPos_real = obj.scenarioObj.sourcesPosition(indRealInActive, :);
            nSrcOrient_real = simulator.vec2rotVec(repmat([0 0 1], [numReal, 1]));
                        
            % Active array loudspeakers
            delays = obj.scenarioObj.delays;
            attenuations = obj.scenarioObj.attenuations;

            numLoudspeakers = obj.scenarioObj.numLoudspeakers;

            if numVirtual == 0
                ldspkrCoef = double.empty(numLoudspeakers, 0);
            else
                ampVirt = obj.amplitude(indVirtual).*obj.virtualVolume(indVirtual);
                phaseVirt = obj.phase(indVirtual);
                freqVirt = obj.frequency(indVirtual);
                delaysVirt = delays(:, indVirtual);
                attenVirt = attenuations(:, indVirtual);
                
                nSrcCoef_virtual = ampVirt.*exp(1i*phaseVirt);
                ldspkrCoef = repmat(nSrcCoef_virtual.', [numLoudspeakers, 1]).*attenVirt.*exp(-1i*2*pi*repmat(freqVirt', [numLoudspeakers, 1]).*delaysVirt);
            end
            
            ldspkrsPos = obj.scenarioObj.loudspeakersPosition;
            ldspkrsOrient = simulator.vec2rotVec(obj.scenarioObj.loudspeakersOrientation);
            
            % Delete the unactive array loudspeakers
            unactiveLdspkrs = false(numLoudspeakers, 1);
            unactiveLdspkrs(obj.channelNumber(indReal)) = true;
            ldspkrCoef = ldspkrCoef(~unactiveLdspkrs, :);
            ldspkrsPos = ldspkrsPos(~unactiveLdspkrs, :);
            ldspkrsOrient = ldspkrsOrient(~unactiveLdspkrs, :);
            numLoudspeakers_new = sum(~unactiveLdspkrs);
                        
            % Unify the coefficients in a single matrix
            coef = zeros(numReal + numLoudspeakers_new, numActiveSour);
            coef(1:numReal, indRealInActive) = diag(nSrcCoef_real);
            coef(numReal+1:end, indVirtualInActive) = ldspkrCoef;
            
            % Set the variables in the simulation object
            obj.simulObj.sourcePositions = [nSrcPos_real; ldspkrsPos];
            obj.simulObj.sourceCoefficients = coef;
            obj.simulObj.sourceOrientations = [nSrcOrient_real; ldspkrsOrient];
            obj.simulObj.radPatFuns = repmat({@(x) simulator.monopoleRadPat(x)}, [numReal + numLoudspeakers_new, 1]);
            obj.simulObj.k = 2*pi*obj.frequency(indActiveSour)/obj.c;

        end
        
        function simulate(obj)
            % Configure simulator object
            obj.createTheoreticalModel();
            
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

