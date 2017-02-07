classdef WFSTool2 < handle
    
    properties
        fig
        player
        reprodPanel
    end
    
    methods
        
        function obj = WFSTool2()
            obj.player = reproductor();
            obj.fig = obj.createFigure();
%             obj.fig = openfig('WFSTool.fig');
%             
%             % Modify original callbacks
%             butTest = findobj(obj.fig.Children, 'Tag', 'but_test');
%             butTest.Callback = @(hObject, eventdata) obj.callback(); %@(hObject, eventdata) but_test_Callback(obj.playObj, Fs);
%             
%             butImportWav = findobj(obj.fig.Children, 'Tag', 'file');
%             table = findobj(obj.fig.Children, 'Tag', 'table');
%             butImportWav.Callback = @(hObject, eventdata) open_Callback(table);
        end
        
        function fig = createFigure(obj)
            fig = figure;
            
            obj.reprodPanel = reproductionPanel(fig, @(action) obj.orderCallback(action));
            
            
        end
             
        function playMusic(obj)
            % Audio file information
            fileName = 'Dangerous Woman.mp3';
            obj.player.audioFileName = fileName;
           
            % Reproduce
            obj.player.executeOrder('play');        
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
                    command.action = 'pause';
                case 'next'
                    % Next track
                    % Get list of tracks
                    list = obj.reprodPanel.list.Data;
                    % Find index of current track
                    [~, index] = ismember(fileName, list);
                    % Increase it by one and check it doesn't fall
                    % out of range
                    N = size(list, 1);
                    index = mod(index, N) + 1;
                    % Get new track
                    activeTrack = list{index};
                case 'previous'
                    % Previous track
                    % Get list of tracks
                    list = obj.reprodPanel.list.Data;
                    % Find index of current track
                    [~, index] = ismember(fileName, list);
                    % Decrease it by one and check it doesn't fall
                    % out of range
                    N = size(list, 1);
                    index = mod(index - 2, N) + 1;
                    % Get new track
                    activeTrack = list{index};                
                case 'doubleClick'
                    % Play other song
                    command.action = 'play';
                    command.fileName = fileName;
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
            end
            
        end
        
    end
    
end

