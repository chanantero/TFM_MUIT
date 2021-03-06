classdef WFSTool2 < handle
    
    properties(SetAccess = private) 
        changed % Variables changed since last update
        virtual
        real
        signalsSpec
        
        fig
        player
        reprodPanel
        scenarioObj
        propPanel
    end
    
    properties(Dependent)
        numSources
    end
       
    % Getters and setters
    methods
        function numSources = get.numSources(obj)
            numSources = numel(obj.signalsSpec);
        end
    end
    
    methods
        
        function obj = WFSTool2()
            fig = figure('Units', 'pixels', 'Position', [0 50 1200 600]);
            obj.fig = fig;
            
            obj.player = reproductor();
            obj.reprodPanel = reproductionPanel(fig, @(action) obj.orderCallback(action));
            obj.scenarioObj = scenario(fig);
            obj.propPanel = propertiesPanel(fig, [0.05 0.1 0.4 0.2]);
            
            obj.changed = struct('virtual', false, 'real', false, 'signalsSpec', false);
                        
            addlistener(obj.player, 'numChannels', 'PostSet', @(~, eventData) obj.changeScenario(eventData.AffectedObject.numChannels(1)));
            addlistener(obj.reprodPanel, 'updatedValues', @(~, evntData) obj.reprodPanelListener(evntData.type));
            addlistener(obj.player, 'playingState', 'PostSet', @(~, eventData) obj.GUIenabling(eventData.AffectedObject.playingState));

            
            obj.signalsSpec = obj.reprodPanel.signals;
            obj.virtual = obj.reprodPanel.virtual;
            obj.real = obj.reprodPanel.real;
            obj.changed.virtual = true;
            obj.changed.real = true;
            obj.changed.signalsSpec = true;
            obj.updateEverything();
        end       
        
        function updateEverything(obj)

            if obj.changed.virtual || obj.changed.real
                rightSize = [obj.numSources, 1];
                virtualRight = all(size(obj.virtual) == rightSize);
                realRight = all(size(obj.real) == rightSize);
                signalsRight = all(size(obj.signalsSpec) == rightSize);
                assert(virtualRight && realRight && signalsRight, 'WFSTool2:updateEverything', 'The signals specifications and the virtual and real flags must have the same size')
                
                obj.updateComMat();
            end
            
            if obj.changed.signalsSpec
                obj.updateSignalProvidersVariables();
            end
        end
        
        function setNumSources(obj, numSources)
            obj.virtual = true(numSources, 1);
            obj.real = true(numSources, 1);
            obj.signalsSpec = cell(numSources, 1);
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
        
    end
    
    methods(Access = private)
        
        function GUIenabling(obj, newPlayingState)
            if newPlayingState == playingStateClass('stopped')
                obj.reprodPanel.enableGUI();
                obj.propPanel.enableGUI();
            else
                obj.reprodPanel.disableGUI();
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
                    
                    obj.changed.virtual = true;
                    obj.changed.real = true;
                    obj.changed.signalsSpec = true;
                case 'virtual'
                    obj.virtual = obj.reprodPanel.virtual;
                    obj.changed.virtual = true;
                case 'real'
                    obj.real = obj.reprodPanel.real;
                    obj.changed.real = true;
            end
            
            obj.updateEverything();
        end
        
        function updateComMat(obj)
            % Update commutation matrix
            comMat = WFSTool2.createCommutationMatrix(obj.virtual, obj.real);
            obj.player.setProps('comMatrix', comMat);
            
            obj.updateSignalProvidersVariables();
            obj.updateProcessorVariables();
            obj.updateGUIConnectionsStuff();
            
            obj.changed.real = false;
            obj.changed.virtual = false;

        end
        
        function updateProcessorVariables(obj)
            % The processors need the getDelayFun and getAttenFun
            % functions, and these are connected to the scenario
            
            % Update scenario based on the virtual sources
            indVirt = find(obj.player.comMatrix(:, 1));
            obj.scenarioObj.setNumSources(numel(indVirt));
            
            % Assign delay and attenuation functions
            for k = 1:numel(indVirt)
                delayFun = @() obj.scenarioObj.delays(:, k);
                obj.player.setProps('getDelayFun', delayFun, [indVirt(k), 1]);
                
                attenFun = @() obj.scenarioObj.attenuations(:, k);
                obj.player.setProps('getAttenFun', attenFun, [indVirt(k), 1]);
            end
            
            % Update the rest of delay and attenuation functions
            [readerIndex, playerIndex] = obj.player.getLinkSubInd();
            % Remove the links of the first column (WFS Array/Virtual sources)
            flag = playerIndex > 1;
            readerIndex = readerIndex(flag);
            playerIndex = playerIndex(flag);
            
            for k = 1:numel(playerIndex)
                numChannels = obj.player.numChannels(playerIndex(k));
                fDelay = @() zeros(numChannels, 1);
                fAtten = @() ones(numChannels, 1);
                obj.player.setProps('getDelayFun', fDelay, [readerIndex(k), playerIndex(k)]);
                obj.player.setProps('getAttenFun', fAtten, [readerIndex(k), playerIndex(k)]);
            end
            
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
        
        function setVolume(obj, volume, indPlayer)
            comMatCoef = obj.player.comMatrixCoef;
            comMatCoef(:, indPlayer) = volume;
            obj.player.setProps('comMatrixCoef', comMatCoef);
        end
        
        function updateSignalProvidersVariables(obj)
            % Signal providers need to have some parameters specified based
            % on the signals specifications
            
            % Assign signals
            for k = 1:obj.numSources
                signalSpec = obj.signalsSpec{k};
                
                % Is a complex number?
                param = regexp(signalSpec, 'A:(?<Amplitude>(\d+\.\d+|\d+)) Ph:(?<Phase>(\d+\.\d+|\d+)) f:(?<Frequency>(\d+\.\d+|\d+))', 'names');
                if isempty(param)
                    % It is not a complex number
                    % Is it a file that exists?
                    a = dir(signalSpec);
                    if numel(a) > 0                        
                        obj.player.setProps('mode', originType('file'), k);
                        obj.player.setProps('audioFileName', signalSpec, k);
                    else
                        % No success. As default, treat it as a tone
                        obj.player.setProps('mode', originType('sinusoidal'), k);
                        obj.player.setProps('amplitude', 0, k);
                        obj.player.setProps('phase', 0, k);
                        obj.player.setProps('frequency', 1, k);
                    end
                else
                    % It is a complex number, set the
                    % properties
                    obj.player.setProps('mode', originType('sinusoidal'), k);
                    obj.player.setProps('amplitude', str2double(param.Amplitude), k);
                    obj.player.setProps('phase', str2double(param.Phase), k);
                    obj.player.setProps('frequency', str2double(param.Frequency), k);
                end
            end
                        
            obj.changed.signalsSpec = false;
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
                case 2
                    sourcePosition = obj.scenarioObj.sourcesPosition;
                    loudspeakersPosition = [-0.1 0 0; 0.1 0 0];
                    loudspeakersOrientation = [1 0 0; -1 0 0];
                    roomPosition = [-2, -2, 4, 4];
                    obj.scenarioObj.setScenario(sourcePosition, loudspeakersPosition, loudspeakersOrientation, roomPosition);
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
                    
                otherwise
                    warning('Wrong number of output channels. There is not possible scenario for that case')
            end
        end
        
    end

    
    methods(Static)
        function comMat = createCommutationMatrix(virtual, real)
            comMat = [virtual, diag(real)];
        end
    end
end

