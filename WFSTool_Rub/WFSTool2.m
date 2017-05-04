classdef WFSTool2 < handle
    
    properties
        fig
        player
        reprodPanel
        scenarioObj
        propPanel
        signalsChanged
        numSources % Number of sources
    end
    
    methods
        
        function obj = WFSTool2()
            fig = figure('Units', 'pixels', 'Position', [0 50 1200 600]);
            obj.fig = fig;
            
            obj.reprodPanel = reproductionPanel(fig, @(action) obj.orderCallback(action));
            obj.scenarioObj = scenario(fig);
            obj.player = reproductor_plus();
            obj.player.set_getDelayFun(@() obj.scenarioObj.delays);
            obj.player.set_getAttenFun(@() obj.scenarioObj.attenuations);
            obj.propPanel = propertiesPanel(obj.player, fig, [0.05 0.1 0.4 0.2], @(x) obj.player.setFrameSize(x), @(x) obj.player.setDevice(x), @(x) obj.player.setDriver(x));
            
            addlistener(obj.player, 'numChannels', 'PostSet', @(~, eventData) obj.changeScenario(eventData.AffectedObject.numChannels));
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
         
    end
    
    methods(Access = private)
        
        function orderCallback(obj, order)
            % Based on the user order and the state of the player, a
            % command for the player is created
            state = obj.player.playingState;
            
            switch order
                case 'play'
                    % Start track again. If signals have changed, assign
                    % the new signals and then reproduce
                    
                    if obj.signalsChanged
                        % Stop the reproductor
                        obj.player.executeOrder('stop');
                        
                        % See if the signals correspond to real
                        % sources, virtual sources or both
                        virtual; % Column vector
                        real; % Column vector
                        
                        % Assign the new commutation matrix
                        comMat = [virtual, diag(double(real))];
                        obj.player.setProps('comMatrix', comMat);
                        
                        % Assign the new signals
                        for k = 1:obj.numSources
                            % Identify if the singal specification is a
                            % file name or a sinusoidal signal
                            signalSpec = obj.reprodPanel.getTrackName(1);
                            
                            % Is a complex number?
                            param = regexp(signalSpec, 'A:(?<Amplitude>(\d+\.\d+|\d+)) Ph:(?<Phase>(\d+\.\d+|\d+)) f:(?<Frequency>(\d+\.\d+|\d+))', 'names');
                            if isempty(param)
                                % It is not a complex number
                                % Is it a file that exists?
                                a = dir(signalSpec);
                                assert(numel(a) > 0, 'The signal specifications must specify an audio file that exists or a complex number');
                                
                                obj.player.setProps('mode', originTypes('file'), k);
                                obj.player.setProps('audioFileName', signalSpec, k);
                            else
                                % It is a complex number, set the
                                % properties
                                obj.player.setProps('mode', originTypes('sinusoidal'), k);
                                obj.player.setProps('amplitude', param.Amplitude, k);
                                obj.player.setProps('phase', param.Phase, k);
                                obj.player.setProps('frequency', param.Frequency, k);
                            end
                            
                        end
                        
                        % Play
                        obj.player.executeOrder('play');
                    else
                        % Just restart
                        obj.player.executeOrder('stop');
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
            
            
            % Send command to the player
            obj.player.executeOrder(command);
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
    
end

