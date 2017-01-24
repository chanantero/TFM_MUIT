classdef reproductionPanel < handle
    
    properties
        panel
        list
    end
    
    methods
        
        function obj = reproductionPanel(parent)
            obj.panel = obj.createPanel(parent);
            obj.list = findobj(obj.panel, 'Tag', 'list');
        end
        
        function panel = createPanel(~, parent)
            panel = uipanel(parent, 'BackgroundColor','white', 'Units', 'normalized', 'Position', [0.05, 0.5, 0.4, 0.4]);
            list = uitable(panel, 'Units', 'Normalized', 'Position', [0.05, 0.2, 0.9, 0.7], ...
                'ColumnName',{'Name'}, 'ColumnFormat', {'char'}, 'Tag', 'list');
            
            addButton = uicontrol(panel, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', [0.05, 0.9, 0.1 0.1]);
            
            y0 = 0.05; h = 0.1; w = 0.1;
            playPos = [0.35, y0, w, h];
            stopPos = [0.2, y0, w, h];
            pausePos = [0.5, y0, w, h];
            prevPos = [0.05, y0, w, h];
            nextPos = [0.65, y0, w, h];
            
            playButton = uicontrol(panel, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', playPos, 'Callback', @(hObject, eventData) obj.computeOrder('play', obj.getState(), obj.getCurrentSong()));
            stopButton = uicontrol(panel, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', stopPos, 'Callback', @(hObject, eventData) obj.computeOrder('stop', obj.getState()));
            pauseButton = uicontrol(panel, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', pausePos, 'Callback', @(hObject, eventData) obj.computeOrder('pause', obj.getState()));
            nextButton = uicontrol(panel, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', nextPos, 'Callback', @(hObject, eventData) obj.computeOrder('next', obj.getState(), obj.getCurrentSong()));
            prevButton = uicontrol(panel, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', prevPos, 'Callback', @(hObject, eventData) obj.computeOrder('previous', obj.getState(), obj.getCurrentSong()));
            
            
            buttons = [addButton, playButton, stopButton, pauseButton, nextButton, prevButton];
            iconNames = {'images/addIcon.png', 'images/playIcon.png', 'images/stopIcon.png', 'images/pauseIcon.jpg', 'images/nextIcon.jpg', 'images/prevIcon.png'};
            
            for k = 1:numel(buttons)
                button = buttons(k);
                imName = iconNames{k};
                
                % Get size in pixels
                button.Units = 'Pixels';
                width = button.Position(3);
                height = button.Position(4);
                % Set image as icon
                icon = imread(imName); % Read image
                icon = imresize(icon, [height, width]); % Resize
                % Convert to binary image
                if size(icon, 3) == 3
                    icon = rgb2gray(icon);
                end
                icon = ind2rgb(icon, [0 0 0; 1 1 1]);
                
                button.CData = icon;
                
                button.Units = 'Normalized';
            end
        end
        
        function state = getState(obj)
        end
        
        function currentSong = getCurrentSong(obj)
        end
        
        function [order, varargout] = computeOrder(obj, action, state, varargin)
            persistent click;
            if isempty(click)
                click = 0;
            end
            switch state
                case playingStateClass('playing')
                    switch action
                        case 'play'
                            % Start song again
                            order = 'play';
                        case 'stop'
                            
                        case 'pause'
                        case 'next'
                        case 'previous'
                        case 'click'
                            if click == 0
                                pause(0.5);
                                click = 0;
                            else click == 1
                                % Double click
                                order = play;
                                song = 
                            end
                                
                    end
                case playingStateClass('stopped')
                    switch action
                        case 'play'
                        case 'stop'
                        case 'pause'
                        case 'next'
                        case 'previous'
                        case 'doubleClick'  
                    end
                case playingStateClass('paused')
                    switch action
                        case 'play'
                        case 'stop'
                        case 'pause'
                        case 'next'
                        case 'previous'
                        case 'doubleClick'
                    end
            end
        end
        
    end
    
end

