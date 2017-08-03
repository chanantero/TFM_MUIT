classdef reproductionPanel_noiseChannel < handle
    
    
    properties(SetAccess = private)
        signals
        virtual
        real
        channelNumber
        numSources
        virtualVolume
        realVolume
    end
    
    properties(SetAccess = private)
        selectedTrack
        panel
        list
        orderCallback
    end
    
    events
        updatedValues
    end
    
    methods
        function obj = reproductionPanel_noiseChannel(parent, position, orderCallback)
            obj.panel = obj.createPanel(parent, position, @(buttonTag) obj.buttonCallback(buttonTag),...
                @(eventData) obj.changeSelection(eventData),...
                @() obj.setTrackCallback(),...
                @(eventData) obj.cellEditCallback(eventData),...
                @(x) obj.numSourcesPopUpMenuCallback(x));
            obj.list = findobj(obj.panel, 'Tag', 'list');
            obj.orderCallback = orderCallback;
            
            obj.setNumSignals(1);
        end
        
        function setSignal(obj, signalSpec, index)
            obj.signals{index} = signalSpec;
            obj.list.Data(index, 1) = {signalSpec};
%             
%             evntData = updatedValuesEvntData('signals');
%             notify(obj, 'updatedValues', evntData);
        end
        
        function setSignals(obj, signalSpecs)
            obj.signals = signalSpecs;
            obj.list.Data(:, 1) = signalSpecs;
            
%             evntData = updatedValuesEvntData('signals');
%             notify(obj, 'updatedValues', evntData);
        end
        
        function setVirtualFlags(obj, flags)
            obj.virtual = flags;
            obj.list.Data(:, 2) = num2cell(flags);
            
%             evntData = updatedValuesEvntData('virtual');
%             notify(obj, 'updatedValues', evntData);
        end
        
        function setRealFlags(obj, flags)
            obj.real = flags;
            obj.list.Data(:, 3) = num2cell(flags);
            
%             evntData = updatedValuesEvntData('real');
%             notify(obj, 'updatedValues', evntData);
        end
        
        function setNumSignals(obj, num)
            obj.list.Data = [cell(num, 1), num2cell(false(num, 2)), num2cell(zeros(num, 1)), ...
                num2cell(ones(num, 1)), num2cell(ones(num, 1))];
            obj.virtual = false(num, 1);
            obj.real = false(num, 1);
            obj.channelNumber = zeros(num, 1);
            obj.virtualVolume = ones(num, 1);
            obj.realVolume = ones(num, 1);
            sign = cell(num, 1);
            for k = 1:num
                sign{k} = '';
            end
            obj.signals = sign;
            
            obj.numSources = num;
            
%             evntData = updatedValuesEvntData('numSources');
%             notify(obj, 'updatedValues', evntData);
        end
        
        function setChannelNumber(obj, channelNumbers)
            obj.channelNumber = channelNumbers;
            obj.list.Data(:, 4) = num2cell(channelNumbers);
        end
        
        function disableGUI(obj)
            obj.list.Enable = 'off';
            numSourcesPopUp = findobj(obj.panel, 'Tag', 'numSources');
            numSourcesPopUp.Enable = 'off';
        end
        
        function enableGUI(obj)
            obj.list.Enable = 'on';
            numSourcesPopUp = findobj(obj.panel, 'Tag', 'numSources');
            numSourcesPopUp.Enable = 'on';
        end
        
    end
    
    methods(Access = private)
        
        function panel = createPanel(~, parent, position, buttonCallback, listSelectionCallback, setTrackCallback, cellEditCallback, setNumSourcesCallback)
            panel = uipanel(parent, 'BackgroundColor','white', 'Units', 'normalized', 'Position', position);
            menu = uicontextmenu;
            list_ = uitable(panel, 'Units', 'Normalized', 'Position', [0.05, 0.2, 0.9, 0.65], ...
                'ColumnName',{'Signal', 'Virtual', 'Real', 'Channel nº', 'Virtual volume', 'Real Volume'}, 'ColumnFormat', {'char', 'logical', 'logical', 'numeric', 'numeric', 'numeric'}, 'Tag', 'list',...
                'ColumnEditable', [true, true, true, true, true, true], 'UIContextMenu', menu,...
                'CellSelectionCallback', @(hObject, eventData) listSelectionCallback(eventData),...
                'CellEditCallback', @(hObject, eventData) cellEditCallback(eventData));
            list_.UserData = struct('paths', []);
            
            uimenu(menu, 'Label', 'Select Audio Track', 'Callback', @(~, ~) setTrackCallback());
            
            y0 = 0.05; h = 0.1; w = 0.1;
            playPos = [0.35, y0, w, h];
            stopPos = [0.2, y0, w, h];
            pausePos = [0.5, y0, w, h];
            
            playButton = uicontrol(panel, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', playPos, 'Callback', @(hObject, eventData) buttonCallback('play'));
            stopButton = uicontrol(panel, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', stopPos, 'Callback', @(hObject, eventData) buttonCallback('stop'));
            pauseButton = uicontrol(panel, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', pausePos, 'Callback', @(hObject, eventData) buttonCallback('pause'));
            
            
            buttons = [playButton, stopButton, pauseButton];
            iconNames = {'images/playIcon.bmp', 'images/stopIcon.bmp', 'images/pauseIcon.jpg'};
            
            for k = 1:numel(buttons)
                button = buttons(k);
                imName = iconNames{k};
                
                % Get size in pixels
                button.Units = 'Pixels';
                width = button.Position(3);
                height = button.Position(4);
                % Set image as icon
                icon = imread(imName); % Read image
                icon = imresize_Rub(icon, [height, width]); % Resize
                % Convert to binary image
                if size(icon, 3) == 3
                    icon = rgb2gray(icon);
                end
                icon = ind2rgb(icon, [0 0 0; 1 1 1]);
                
                button.CData = icon;
                
                button.Units = 'Normalized';
            end
            
            % Number of sources uicontrol
            uicontrol(panel, 'Style', 'text', 'String', 'Nº sources', 'Units', 'normalized',...
                'Position', [0 0.9 0.2 0.1]);
            uicontrol(panel, 'Style', 'popupmenu', 'Units', 'normalized', 'Position', [0.2 0.9 0.1 0.1],...
                'Tag', 'numSources', 'String', {'1', '2', '3'}, 'Callback', @(hObject, ~, ~) setNumSourcesCallback(str2double(hObject.String{hObject.Value})));
        end
        
        function buttonCallback(obj, buttonTag)
            % Specify what the order is
            switch buttonTag
                case 'play'
                    order = 'play';
                case 'stop'
                    order = 'stop';
                case 'pause'
                    order = 'pause';
            end
            
            % Send order to the callback
            obj.orderCallback(order);
        end
        
        function numSourcesPopUpMenuCallback(obj, numSources)
            obj.setNumSignals(numSources);
        end
        
        function changeSelection(obj, eventdata)
            if ~isempty(eventdata.Indices)
                obj.selectedTrack = eventdata.Indices(1);
            else
                obj.selectedTrack = [];
            end
        end
        
        function setTrackCallback(obj)
            musicDirectory = 'C:\Users\Rubén\Music\Varias\';
            [FileName, PathName, ~] = uigetfile({'*.mp3', 'MP3 File'; '*.wav', 'WAV File'}, 'Select audio track', 'MultiSelect', 'on', musicDirectory);
            
            
            if FileName ~= 0
                obj.setSignal([PathName, FileName], obj.selectedTrack);
            end
        end
        
        function cellEditCallback(obj, eventData)
            
            row = eventData.Indices(1);
            column = eventData.Indices(2);
            
            newSignal = eventData.EditData;
            
            if column == 1
                % If the name of the track has been changed, check:
                % Is a complex number?
                param = regexp(newSignal, 'A:(?<Amplitude>(\d+\.\d+|\d+)) Ph:(?<Phase>(\d+\.\d+|\d+)) f:(?<Frequency>(\d+\.\d+|\d+))', 'names');
                a = dir(newSignal);
                complexFlag = ~isempty(param);
                fileFlag = numel(a) > 0;
                if complexFlag == false && fileFlag == false
                    % Wrong format
                    % Restart
                    obj.list.Data{row, 1} = 'A:0 Ph:0 f:1';
                end
                obj.signals{row} = newSignal;
                
                evntData = updatedValuesEvntData('signals');
                notify(obj, 'updatedValues', evntData);
            elseif column == 2
                % The virtual flags have been changed
                obj.virtual = cell2mat(obj.list.Data(:, 2));
                
                evntData = updatedValuesEvntData('virtual');
                notify(obj, 'updatedValues', evntData);
            elseif column == 3
                % The real flags have been changed
                obj.real = cell2mat(obj.list.Data(:, 3));
                
                evntData = updatedValuesEvntData('real');
                notify(obj, 'updatedValues', evntData);
            elseif column == 4
                % The number of the channels for the sources have changed
                candidate = cell2mat(obj.list.Data(:, 4));
                if sum(obj.real) > numel(unique(candidate(obj.real)));
                    % Repeated channel
                    obj.list.Data(row, column) = eventData.PreviousData;
                else
                    obj.channelNumber = candidate;
                end
                
                evntData = updatedValuesEvntData('channelNumber');
                notify(obj, 'updatedValues', evntData);
            elseif column == 5
                obj.virtualVolume = cell2mat(obj.list.Data(:, 5));
                
                evntData = updatedValuesEvntData('virtualVolume');
                notify(obj, 'updatedValues', evntData);
            elseif column == 6
                obj.realVolume = cell2mat(obj.list.Data(:, 6));
                
                evntData = updatedValuesEvntData('realVolume');
                notify(obj, 'updatedValues', evntData);
            end
            
        end
        
    end
    
    
end

