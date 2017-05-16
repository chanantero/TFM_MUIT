classdef WFSTool2 < handle
    
    properties(SetAccess = private) 
        changed % Variables changed since last update
        numSources
        virtual
        real
        signalsSpec
        
        fig
        player
        reprodPanel
        scenarioObj
        propPanel
    end
        
    methods
        
        function obj = WFSTool2()
            fig = figure('Units', 'pixels', 'Position', [0 50 1200 600]);
            obj.fig = fig;
            
            obj.player = reproductor();
            obj.reprodPanel = reproductionPanel(fig, @(action) obj.orderCallback(action));
            obj.scenarioObj = scenario(fig);
            obj.propPanel = propertiesPanel(obj.player, fig, [0.05 0.1 0.4 0.2]);
            
            obj.changed = struct('numSources', false, 'virtual', false, 'real', false, 'signalsSpec', false);
            
            obj.setNumSources(1);
            obj.updateEverything();
                        
            addlistener(obj.player, 'numChannels', 'PostSet', @(~, eventData) obj.changeScenario(eventData.AffectedObject.numChannels(1)));
            addlistener(obj.reprodPanel, 'numSources', 'PostSet', @(~, ~) obj.reprodPanelListener('numSources'));
            addlistener(obj.reprodPanel, 'virtual', 'PostSet', @(~, ~) obj.reprodPanelListener('virtual'));
            addlistener(obj.reprodPanel, 'real', 'PostSet', @(~, ~) obj.reprodPanelListener('real'));
            addlistener(obj.reprodPanel, 'signals', 'PostSet', @(~, ~) obj.reprodPanelListener('signals'));            
        end       
        
        function noRealTimeReproduction(obj, t, positions)
            % Calculate the delays and attenuations for each position
            % Use scenario object
            numChann = obj.player.numChannels;
            delay = zeros(numel(t), numChann);
            attenuation = zeros(numel(t), numChann);
            for k = 1:numel(t)
                % Set scenario
                obj.scenarioObj.setSourcePosition(positions(k, :));
                
                % Get delay and attenuation of the different channels
                delay(k, :) = obj.scenarioObj.getDelays();
                attenuation(k, :) = obj.scenarioObj.getAttenuations();
            end
            
            obj.player.reproduceNoRealTime(t, delay, attenuation);
        end
        
        function updateEverything(obj)
            if obj.changed.numSources
                obj.updateNumSourcesStuff();
            else
                if obj.changed.virtual || obj.changed.real
                    obj.updateComMat();
                end
                
                if obj.changed.signalsSpec
                    obj.updateSignalProvidersVariables();
                end
            end
        end
        
        function setNumSources(obj, value)
            obj.numSources = value;
            obj.changed.numSources = true;
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
        
        function reprodPanelListener(obj, name)
            switch name
                case 'signals'
                    obj.signalsSpec = obj.reprodPanel.signals;
                    obj.changed.signalsSpec = true;
                case 'numSources'
                    obj.numSources = obj.reprodPanel.numSources;
%                     obj.virtual = obj.reprodPanel.virtual;
%                     obj.real = obj.reprodPanel.real;
%                     obj.signalsSpec = obj.reprodPanel.signals;
                    
%                     obj.changed.virtual = true;
%                     obj.changed.real = true;
%                     obj.changed.signalsSpec = true;
                    obj.changed.numSources = true;
                case 'virtual'
                    obj.virtual = obj.reprodPanel.virtual;
                    obj.changed.virtual = true;
                case 'real'
                    obj.real = obj.reprodPanel.real;
                    obj.changed.real = true;
            end
            
            obj.updateEverything();
        end
        
        function updateNumSourcesStuff(obj)
            N = obj.numSources;
            obj.virtual = true(N, 1);
            obj.real = true(N, 1);
            obj.signalsSpec = cell(N, 1);
            for k = 1:N
                obj.signalsSpec{k} = '';
            end
            
            obj.updateComMat();
            obj.updateSignalProvidersVariables();     
            
            obj.changed.numSources = false;
        end
        
        function updateComMat(obj)
            % Update commutation matrix
            comMat = WFSTool2.createCommutationMatrix(obj.virtual, obj.real);
            obj.player.setProps('comMatrix', comMat);
            
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
                delayFun = @() obj.scenarioObj.delays(:, indVirt(k));
                obj.player.setProps('getDelayFun', delayFun, [indVirt(k), 1]);
                
                attenFun = @() obj.scenarioObj.attenuations(:, indVirt(k));
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
            indReal = [1, find(any(obj.player.comMatrix(:, 2:end), 1))+1];
            setFrameDurationFunc = @(frameDuration) obj.player.setProps('frameDuration', frameDuration);
            setDeviceFunc = @(device, index) obj.player.setProps('device', device, indReal(index));
            setDriverFunc = @(driver, index) obj.player.setProps('driver', driver, indReal(index));
            getAvailableDevicesFunc = @(index) obj.player.player{indReal(index)}.getAvailableDevices();
            labels = cell(numel(indReal), 1);
            labels{1} = 'WFS Array';
            for k = 2:numel(indReal)
                labels{k} = sprintf('Device %d', indReal(k)-1);
            end

            obj.propPanel.setFunctions(setFrameDurationFunc, setDeviceFunc, setDriverFunc, getAvailableDevicesFunc, labels);
            
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
                        warning('WFSTool2:updateSignalProvidersVariables', 'The signal specifications must specify an audio file that exists or a complex number');
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
            switch numChannels
                case 0;
                    sourcePosition = [0 0 0];
                    loudspeakersPosition = double.empty(0,3);
                    loudspeakersOrientation = double.empty(0,3);
                    roomPosition = [0 0 1 1];
                    obj.scenarioObj.setScenario(sourcePosition, loudspeakersPosition, loudspeakersOrientation, roomPosition);
                case 2;
                    sourcePosition = [0 1 0];
                    loudspeakersPosition = [-0.1 0 0; 0.1 0 0];
                    loudspeakersOrientation = [1 0 0; -1 0 0];
                    roomPosition = [-2, -2, 4, 4];
                    obj.scenarioObj.setScenario(sourcePosition, loudspeakersPosition, loudspeakersOrientation, roomPosition);
                case 96;                    
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

