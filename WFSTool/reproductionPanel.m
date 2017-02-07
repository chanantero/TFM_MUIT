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
            obj.panel = obj.createPanel(parent, @(buttonTag) obj.buttonCallback(buttonTag),...
                @() obj.listCallback(),...
                @(eventData) obj.changeSelection(eventData),...
                @() obj.addTrackCallback());
            obj.list = findobj(obj.panel, 'Tag', 'list');
            obj.orderCallback = orderCallback;
        end
        
        function panel = createPanel(~, parent, buttonCallback, listButtonDownCallback, listSelectionCallback, addButtonCallback)
            panel = uipanel(parent, 'BackgroundColor','white', 'Units', 'normalized', 'Position', [0.05, 0.5, 0.4, 0.4]);
            list = uitable(panel, 'Units', 'Normalized', 'Position', [0.05, 0.2, 0.9, 0.7], ...
                'ColumnName',{'Name'}, 'ColumnFormat', {'char'}, 'Tag', 'list',...
                'ButtonDownFcn', @(hObject, eventData) listButtonDownCallback(),...
                'CellSelectionCallback', @(hObject, eventData) listSelectionCallback(eventData));
            list.UserData = struct('paths', []);
            
            addButton = uicontrol(panel, 'Style', 'pushbutton', 'Tag', 'addButton',...
                'Units', 'normalized', 'Position', [0.05, 0.9, 0.1 0.1], 'Callback', @(hObject, eventData) addButtonCallback());
            
            y0 = 0.05; h = 0.1; w = 0.1;
            playPos = [0.35, y0, w, h];
            stopPos = [0.2, y0, w, h];
            pausePos = [0.5, y0, w, h];
            prevPos = [0.05, y0, w, h];
            nextPos = [0.65, y0, w, h];
            
            playButton = uicontrol(panel, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', playPos, 'Callback', @(hObject, eventData) buttonCallback('play'));
            stopButton = uicontrol(panel, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', stopPos, 'Callback', @(hObject, eventData) buttonCallback('stop'));
            pauseButton = uicontrol(panel, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', pausePos, 'Callback', @(hObject, eventData) buttonCallback('pause'));
            nextButton = uicontrol(panel, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', nextPos, 'Callback', @(hObject, eventData) buttonCallback('next'));
            prevButton = uicontrol(panel, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', prevPos, 'Callback', @(hObject, eventData) buttonCallback('previous'));
            
            
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
        
        function buttonCallback(obj, buttonTag)
            % Create order structure
            switch buttonTag
                case 'play'
                    order.action = 'play';
                case 'stop'
                    order.action = 'stop';
                case 'pause'
                    order.action = 'pause';
                case 'next'
                    order.action = 'next';
                    order.fileName = 'activeTrack';
                case 'previous'
                    order.action = 'previous';
                    order.fileName = 'activeTrack';
            end
            
            % Send order structure to the callback
            obj.orderCallback(order);
        end
        
        function listCallback(obj)
            persistent click lastClickedSong;
            if isempty(click)
                click = 0;
            end
            if isempty(lastClickedSong)
                lastClickedSong = '';
            end
            
            if click == 0
                % Single click.
                click = 1;
                obj.list.Enable = 'innactive';
                pause(0.5);
                click = 0;
                obj.list.Enable = 'on';
            else
                % Double click                
                % Create order structure
                order.action = 'play';
                order.fileName = obj.selectedTrack;
                click = 0;
                obj.list.Enable = 'on';
               
                % Send order structure to the callback
                obj.orderCallback(order);
            end       
            
        end
        
        function changeSelection(obj, eventdata)
            obj.selectedTrack = obj.list.Data{eventdata.Indices(1)};
        end
        
        function setSelectedTrack(obj, trackName)
            obj.selectedTrack = trackName;
        end
        
        function setActiveTrack(obj, trackName)
            obj.activeTrack = trackName;
        end
        
        function addTrackCallback(obj)
            [FileName, PathName, ~] = uigetfile({'*.wav', 'WAV File'; '*.mp3', 'MP3 File'}, 'Select audio track', 'MultiSelect', 'on');
            
            if ~iscell(FileName)
                FileName = {FileName};
                PathName = {PathName};
            end
            
            for k = 1:numel(FileName)
                obj.list.Data = [obj.list.Data; FileName(k)];
                obj.list.UserData.paths = [obj.list.UserData.paths; PathName(k)];
            end
        end
        
    end
    
end

