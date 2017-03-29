classdef WFSTool2 < handle
    
    properties
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
            
            obj.reprodPanel = reproductionPanel(fig, @(action) obj.orderCallback(action));
            obj.scenarioObj = scenario(fig);
            obj.player = reproductor();
            obj.player.set_getDelayFun(@() obj.scenarioObj.delays);
            obj.player.set_getAttenFun(@() obj.scenarioObj.attenuations);
            obj.propPanel = propertiesPanel(obj.player, fig, [0.05 0.1 0.4 0.2], @(x) obj.player.setFrameSize(x), @(x) obj.player.setDevice(x), @(x) obj.player.setDriver(x));
            
            addlistener(obj.player, 'numChannels', 'PostSet', @(~, eventData) obj.changeScenario(eventData.AffectedObject.numChannels));
        end       
        
    end
    
    methods(Access = private)
        function command = createCommand(obj, order, state)
            
            p = inputParser;
            addParameter(p, 'action', '')
            addParameter(p, 'fileName', '')
            parse(p, order)
            
            action = p.Results.action;
            fileName = p.Results.fileName;
            
            switch action
                case 'play'
                    % Start track again
                    command.action = 'play';
                case 'stop'
                    % Stop
                    command.action = 'stop';
                case 'pause'
                    % Pause
                    switch state
                        case playingStateClass('playing')
                            command.action = 'pause';
                        case playingStateClass('stopped')
                            % Do nothing
                            command = [];
                        case playingStateClass('paused')
                            % Resume
                            command.action = 'resume';
                    end
                case 'next'
                    % Next track
%                     % Get list of tracks
%                     list = obj.reprodPanel.getTrackNames();
%                     % Find index of current track
%                     [~, index] = ismember(fileName, list);
%                     % Increase it by one and check it doesn't fall
%                     % out of range
%                     N = size(list, 1);
%                     index = mod(index, N) + 1;
%                     % Get new track
%                     activeTrack = list{index};
                    % Based in the active track registered in the
                    % reproduction panel object
                    activeTrack = obj.reprodPanel.getNextTrackName();
                case 'previous'
                    % Previous track
                    activeTrack = obj.reprodPanel.getPreviousTrackName();       
                case 'doubleClick'
                    % Play other song
                    command.action = 'play';
                    command.fileName = fileName;
                    % Update panel
                    obj.reprodPanel.setActiveTrack(fileName);
            end
            
            if ismember(action, {'next', 'previous'})
                switch state
                    case playingStateClass('playing')
                        command.action = 'play';
                    case playingStateClass('stopped')
                        command.action = 'assignTrack';
                    case playingStateClass('paused')
                        command.action = 'assignTrack';
                end
                command.fileName = activeTrack;
                % Update panel
                obj.reprodPanel.setActiveTrack(activeTrack);
            end
            
        end
        
        function orderCallback(obj, order)
            % order is a structure with the information of the order the
            % user has specified when interacting with the reproduction
            % panel
            
            % Based on the user order and the state of the player, a
            % command for the player is created
            state = obj.player.playingState;
            command = obj.createCommand(order, state);
            
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

