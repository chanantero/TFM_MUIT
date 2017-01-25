classdef reproductionPanel < handle
    
    properties
        panel
        list
        orderCallback
    end
    
    properties(SetAccess = private)
        selectedTrack
        activeTrack
    end
    
    methods
        
        function obj = reproductionPanel(parent, orderCallback)
            obj.panel = obj.createPanel(parent, buttonCallback);
            obj.list = findobj(obj.panel, 'Tag', 'list');
            obj.orderCallback = orderCallback;
        end
        
        function panel = createPanel(~, parent, buttonCallback, listCallback)
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
                'Units', 'normalized', 'Position', playPos, 'Callback', @(hObject, eventData) buttonCallback(struct('action', 'play')));
            stopButton = uicontrol(panel, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', stopPos, 'Callback', @(hObject, eventData) buttonCallback(struct('action', 'stop')));
            pauseButton = uicontrol(panel, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', pausePos, 'Callback', @(hObject, eventData) buttonCallback(struct('action', 'pause')));
            nextButton = uicontrol(panel, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', nextPos, 'Callback', @(hObject, eventData) buttonCallback(struct('action', 'next')));
            prevButton = uicontrol(panel, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', prevPos, 'Callback', @(hObject, eventData) buttonCallback(struct('action', 'previous')));
            
            
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
        
        function setSelectedTrack(~, trackName)
            
            obj.selectedTrack = trackName;
        end
        
        function setActiveTrack(~, trackName)
            
        end
        
    end
    
end

