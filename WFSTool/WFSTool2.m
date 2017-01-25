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
        
        function orderCallback(obj, action)
            
        end
       
        function command = createCommand(obj, order, state, activeTrack)
            persistent click lastClickedSong;
            if isempty(click)
                click = 0;
            end
            if isempty(lastClickedSong)
                lastClickedSong = '';
            end
            
            p = inputParser;
            addParameter(p, 'action', '')
            addParameter(p, 'fileName', '')
            parse(p, order)
            
            action = p.Results.action;
            fileName = p.Results.fileName;
            
            switch state
                case playingStateClass('playing')
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
                            [~, index] = ismember(activeTrack, list);
                            % Increase it by one and check it doesn't fall
                            % out of range
                            N = size(list, 1);
                            index = mod(index, N) + 1;
                            % Get new track
                            activeTrack = list{index};
                            % Create order structure
                            command.action = 'play';
                            command.fileName = activeTrack;
                        case 'previous'
                            % Previous track
                            % Get list of tracks
                            list = obj.reprodPanel.list.Data;
                            % Find index of current track
                            [~, index] = ismember(activeTrack, list);
                            % Decrease it by one and check it doesn't fall
                            % out of range
                            N = size(list, 1);
                            index = mod(index - 2, N) + 1;
                            % Get new track
                            activeTrack = list{index};
                            % Create order structure
                            command.action = 'play';
                            command.fileName = activeTrack;
                        case 'doubleClick'
                            % Play other song
                            command.action = 'play';
                            command.fileName = fileName;
%                             if click == 0
%                                 click = 1;
%                                 pause(0.5);
%                                 click = 0;
%                             else click == 1
%                                 % Double click
%                                 % Get The file
%                                 order.action = play;
%                                 order.fileName = stop;
%                             end
                                
                    end
                case playingStateClass('stopped')
                    switch action
                        case 'play'
                            % Play active track
                            command.action = 'play';
                        case 'stop'
                            % Stop (do nothing)
                            command.action = 'stop';
                        case 'pause'
                            % Pause (do nothing)
                            command.action = 'pause';
                        case 'next'
                            % Next track
                            % Get list of tracks
                            list = obj.reprodPanel.list.Data;
                            % Find index of current track
                            [~, index] = ismember(activeTrack, list);
                            % Increase it by one and check it doesn't fall
                            % out of range
                            N = size(list, 1);
                            index = mod(index, N) + 1;
                            % Get new track
                            activeTrack = list{index};
                            % Create order structure                           
                            command.action = 'assignTrack';
                            command.fileName = activeTrack;                         
                        case 'previous'
                            % Previous track
                            % Get list of tracks
                            list = obj.reprodPanel.list.Data;
                            % Find index of current track
                            [~, index] = ismember(activeTrack, list);
                            % Decrease it by one and check it doesn't fall
                            % out of range
                            N = size(list, 1);
                            index = mod(index - 2, N) + 1;
                            % Get new track
                            activeTrack = list{index};
                            % Create order structure
                            command.action = 'assignTrack';
                            command.fileName = activeTrack;
                        case 'doubleClick'
                            % Play other song
                            command.action = 'play';
                            command.fileName = fileName;
                    end
                case playingStateClass('paused')
                    switch action
                        case 'play'
                        case 'stop'
                        case 'pause'
                        case 'next'
                        case 'previous'
                        case 'doubleClick'
                            % Play other song
                            command.action = 'play';
                            command.fileName = fileName;
                    end
            end
        end
        
    end
    
end

