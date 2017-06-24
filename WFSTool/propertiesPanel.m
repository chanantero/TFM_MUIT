classdef propertiesPanel < handle
    
    properties
        panel
        menuLabels
        frameDurationMenu
        fieldLabels
        volumeControls
        deviceMenus
        driverMenus
        
        readerMenuLabels
        readerDriverMenus
        readerDeviceMenus
    end
    
    methods
        function obj = propertiesPanel(fig, position, setFrameDurationFunc, setVolumeFunc, setDeviceFunc, setDriverFunc, getAvailableDevicesFunc, labels)            
            obj.panel = uipanel(fig, 'Units', 'normalized', 'Position', position);
            
            if nargin == 8
                obj.setFunctionsWriters(setFrameDurationFunc, setVolumeFunc, setDeviceFunc, setDriverFunc, getAvailableDevicesFunc, labels);
            end
            
        end
        
        function setFunctionsReaders(obj, setDeviceFunc, setDriverFunc, getAvailableDevicesFunc, labels)
            % Delete previous menus
            delete(obj.readerDeviceMenus);
            delete(obj.readerDriverMenus);
            delete(obj.readerMenuLabels);
                
            % Create new ones
            N = numel(labels);
            obj.readerDeviceMenus = gobjects(N, 1);
            obj.readerDriverMenus = gobjects(N, 1);
            obj.readerMenuLabels = gobjects(N, 1);
                        
            % Positions
            gridPos = [0.15 0 0.65 0.3];
            numCol = 3; % Number of fields
            numRow = N;
            xRelMargin = 0.1; % 10%
            yRelMargin = 0.1; % 10%
            
            squareXDim = gridPos(3)/numCol;
            squareYDim = gridPos(4)/max(numRow, 3); 
            squaresXPos = gridPos(1) + squareXDim*(0:numCol-1);
            squaresYPos = gridPos(2) + gridPos(4) - squareYDim*(1:numRow);
            % Apply margins
            xMargin = xRelMargin*squareXDim;
            yMargin = yRelMargin*squareYDim;    
            squaresXPos = squaresXPos + xMargin;
            squaresYPos = squaresYPos + yMargin;
            squareXDim = squareXDim - 2*xMargin;
            squareYDim = squareYDim - 2*yMargin;
            
            xDim_label = 0.15;
            xPos_label = squaresXPos(1) - xDim_label;
            yPos = squaresYPos;
            yDim = squareYDim;
            xPos_deviceMenu = squaresXPos(2);
            xDim_deviceMenu = squareXDim;
            xPos_driverMenu = squaresXPos(3);
            xDim_driverMenu = squareXDim;
                       
            for k = 1:N           
                obj.readerDeviceMenus(k) = uicontrol(obj.panel, 'Style', 'popupmenu', 'String', getAudioDevices(audioDeviceWriter()), 'Tag', 'device', 'Units', 'normalized', 'Position', [xPos_deviceMenu, yPos(k), xDim_deviceMenu, yDim],...
                    'Callback', @(hObject, ~) setDeviceFunc(hObject.String{hObject.Value}, k));
                
                obj.readerDriverMenus(k) = uicontrol(obj.panel, 'Style', 'popupmenu', 'String', {'DirectSound', 'ASIO'}, 'Tag', 'driver', 'Units', 'normalized', 'Position', [xPos_driverMenu, yPos(k), xDim_driverMenu, yDim],...
                    'Callback', @(hObject, ~) driverMenuCallback(hObject.String{hObject.Value}, k));
                
                obj.readerMenuLabels(k) = uicontrol(obj.panel, 'Style', 'text', 'String', labels{k}, 'Tag', 'label', 'Units', 'normalized', 'Position', [xPos_label, yPos(k), xDim_label, yDim]);
            end
            
            % Execute callbacks for the first time
            
            for k = 1:N
                driverOptions = obj.readerDriverMenus(k).String;
                driver = driverOptions{obj.readerDriverMenus(k).Value};
                driverMenuCallback(driver, k);
                
                deviceOptions = obj.readerDeviceMenus(k).String;
                device = deviceOptions{obj.readerDeviceMenus(k).Value};
                setDeviceFunc(device, k);      
            end
            
            
            function driverMenuCallback(driver, ind)
                % Update deviceMenu
                setDriverFunc(driver, ind);
                
                % Get devices for the new driver
                obj.readerDeviceMenus(ind).String = getAvailableDevicesFunc(ind);
                obj.readerDeviceMenus(ind).Value = 1;
            end
            
        end
        
        function setFunctionsWriters(obj, setFrameDurationFunc, setVolumeFunc, setDeviceFunc, setDriverFunc, getAvailableDevicesFunc, labels)
            % Delete previous menus
            delete(obj.frameDurationMenu);
            delete(obj.volumeControls);
            delete(obj.deviceMenus);
            delete(obj.driverMenus);
            delete(obj.menuLabels);
            delete(obj.fieldLabels);
                
            % Create new ones
            N = numel(labels);
            obj.frameDurationMenu = gobjects(1, 1);
            obj.volumeControls = gobjects(N, 1);
            obj.deviceMenus = gobjects(N, 1);
            obj.driverMenus = gobjects(N, 1);
            obj.menuLabels = gobjects(N, 1);
            obj.fieldLabels = gobjects(3, 1);
                        
            % Positions
            gridPos = [0.15 0.3 0.65 0.3];
            numCol = 3; % Number of fields
            numRow = N;
            xRelMargin = 0.1; % 10%
            yRelMargin = 0.1; % 10%
            
            squareXDim = gridPos(3)/numCol;
            squareYDim = gridPos(4)/max(numRow, 3); 
            squaresXPos = gridPos(1) + squareXDim*(0:numCol-1);
            squaresYPos = gridPos(2) + gridPos(4) - squareYDim*(1:numRow);
            % Apply margins
            xMargin = xRelMargin*squareXDim;
            yMargin = yRelMargin*squareYDim;    
            squaresXPos = squaresXPos + xMargin;
            squaresYPos = squaresYPos + yMargin;
            squareXDim = squareXDim - 2*xMargin;
            squareYDim = squareYDim - 2*yMargin;
            
            yPos_frameDurationMenu = 0.92;
            yDim_frameDurationMenu = 0.1;
            xDim_label = 0.15;
            xPos_label = squaresXPos(1) - xDim_label;
            yPos_FieldLabels = gridPos(2) + gridPos(4);
            yDim_FieldLabels = 0.1;
            yPos = squaresYPos;
            yDim = squareYDim;
            xPos_Volume = squaresXPos(1);
            xDim_Volume = squareXDim;
            xPos_deviceMenu = squaresXPos(2);
            xDim_deviceMenu = squareXDim;
            xPos_driverMenu = squaresXPos(3);
            xDim_driverMenu = squareXDim;
            
            
            obj.fieldLabels(1) = uicontrol(obj.panel, 'Style', 'text', 'String', 'Volume', 'Tag', 'FieldLabel', 'Units', 'normalized', 'Position', [xPos_Volume, yPos_FieldLabels, xDim_Volume, yDim_FieldLabels]);
            obj.fieldLabels(2) = uicontrol(obj.panel, 'Style', 'text', 'String', 'Device', 'Tag', 'FieldLabel', 'Units', 'normalized', 'Position', [xPos_deviceMenu, yPos_FieldLabels, xDim_deviceMenu, yDim_FieldLabels]);
            obj.fieldLabels(3) = uicontrol(obj.panel, 'Style', 'text', 'String', 'Driver', 'Tag', 'FieldLabel', 'Units', 'normalized', 'Position', [xPos_driverMenu, yPos_FieldLabels, xDim_driverMenu, yDim_FieldLabels]);
            
          
            uicontrol(obj.panel, 'Style', 'text', 'String', 'Frame Duration (s)', 'Tag', 'frameDurationLabel', 'Units', 'normalized', 'Position', [0.05, yPos_frameDurationMenu, 0.25, yDim_frameDurationMenu]);
            obj.frameDurationMenu = uicontrol(obj.panel, 'Style', 'popupmenu', 'String', {'0.5', '1', '2'}, 'Tag', 'frameDuration',...
                    'Units', 'normalized', 'Position', [0.3, yPos_frameDurationMenu, 0.15, yDim_frameDurationMenu], 'Callback', @(hObject, ~) setFrameDurationFunc(str2double(hObject.String{hObject.Value})));
             
            for k = 1:N     
                obj.volumeControls(k) = uicontrol(obj.panel, 'Style', 'Edit', 'String', '1', 'Tag', 'Volume', 'Units', 'normalized', 'Position', [xPos_Volume, yPos(k), xDim_Volume, yDim],...
                    'Callback', @(hObject, ~) setVolumeFunc(str2double(hObject.String), k));
                
                obj.deviceMenus(k) = uicontrol(obj.panel, 'Style', 'popupmenu', 'String', getAudioDevices(audioDeviceWriter()), 'Tag', 'device', 'Units', 'normalized', 'Position', [xPos_deviceMenu, yPos(k), xDim_deviceMenu, yDim],...
                    'Callback', @(hObject, ~) setDeviceFunc(hObject.String{hObject.Value}, k));
                
                obj.driverMenus(k) = uicontrol(obj.panel, 'Style', 'popupmenu', 'String', {'DirectSound', 'ASIO'}, 'Tag', 'driver', 'Units', 'normalized', 'Position', [xPos_driverMenu, yPos(k), xDim_driverMenu, yDim],...
                    'Callback', @(hObject, ~) driverMenuCallback(hObject.String{hObject.Value}, k));
                
                obj.menuLabels(k) = uicontrol(obj.panel, 'Style', 'text', 'String', labels{k}, 'Tag', 'label', 'Units', 'normalized', 'Position', [xPos_label, yPos(k), xDim_label, yDim]);
            end
            
            % Execute callbacks for the first time
            frameDuration = str2double(obj.frameDurationMenu.String(obj.frameDurationMenu.Value));
            setFrameDurationFunc(frameDuration);
            
            for k = 1:N
                driverOptions = obj.driverMenus(k).String;
                driver = driverOptions{obj.driverMenus(k).Value};
                driverMenuCallback(driver, k);
                
                deviceOptions = obj.deviceMenus(k).String;
                device = deviceOptions{obj.deviceMenus(k).Value};
                setDeviceFunc(device, k);      
                
                volume = str2double(obj.volumeControls(k).String);
                setVolumeFunc(volume, k);
            end
            
            
            function driverMenuCallback(driver, ind)
                % Update deviceMenu
                setDriverFunc(driver, ind);
                
                % Get devices for the new driver
                obj.deviceMenus(ind).String = getAvailableDevicesFunc(ind);
                obj.deviceMenus(ind).Value = 1;
            end
            
        end
        
        function enableGUI(obj)
            for k = 1:numel(obj.driverMenus)
                obj.frameDurationMenu.Visible = 'on';
%                 obj.volumeControls(k).Visible = 'on';
                obj.deviceMenus(k).Visible = 'on';
                obj.driverMenus(k).Visible = 'on';
%                 obj.menuLabels(k).Visible = 'on';
            end
            
            for k = 1:numel(obj.driverMenus)
                obj.readerDeviceMenus(k).Visible = 'on';
                obj.readerDriverMenus(k).Visible = 'on';
            end
        end
        
        function disableGUI(obj)
            for k = 1:numel(obj.driverMenus)
                obj.frameDurationMenu.Visible = 'off';
%                 obj.volumeControls(k).Visible = 'off';
                obj.deviceMenus(k).Visible = 'off';
                obj.driverMenus(k).Visible = 'off';
%                 obj.menuLabels(k).Visible = 'off';
            end
            
            for k = 1:numel(obj.driverMenus)
                obj.readerDeviceMenus(k).Visible = 'off';
                obj.readerDriverMenus(k).Visible = 'off';
            end
        end
       
    end
    
    
end


