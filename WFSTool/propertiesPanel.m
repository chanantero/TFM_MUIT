classdef propertiesPanel < handle
    
    properties
        panel
        frameSizeMenu
        deviceMenu
        driverMenu
        reprodObj
    end
    
    methods
        function obj = propertiesPanel(reprodObj, fig, position, setFrameSize, setDevice, setDriver)
            obj.panel = obj.createGraphics(fig, position, setFrameSize, setDevice, @driverMenuCallback);
                    
            function driverMenuCallback(driver)
                % Update deviceMenu
                setDriver(driver);
                
                % Get devices for the new driver
                obj.deviceMenu.String = reprodObj.getWritingDevices();
                obj.deviceMenu.Value = 1;
            end
            
            addlistener(reprodObj, 'playingState', 'PostSet', @(~, eventData) obj.playingStateSetCallback(eventData.AffectedObject.playingState));
        end
    end
    
    methods(Access = private)
        function panel = createGraphics(obj, parent, position, setFrameSize, setDevice, setDriver)
            %             position = [0.1, 0.1, 0.8, 0.8];
            panel = uipanel(parent, 'Units', 'normalized', 'Position', position);
            
            obj.frameSizeMenu = uicontrol(panel, 'Style', 'popupmenu', 'String', {'1024', '4096', '16384'}, 'Tag', 'frameSize',...
                'Units', 'normalized', 'Position', [0.1, 0.1, 0.2, 0.8], 'Callback', @(hObject, ~) setFrameSize(str2double(hObject.String{hObject.Value})));
            
            obj.deviceMenu = uicontrol(panel, 'Style', 'popupmenu', 'String', getAudioDevices(audioDeviceWriter()), 'Tag', 'device', 'Units', 'normalized', 'Position', [0.4, 0.1, 0.2, 0.8],...
                'Callback', @(hObject, ~) setDevice(hObject.String{hObject.Value}));
            
            obj.driverMenu = uicontrol(panel, 'Style', 'popupmenu', 'String', {'DirectSound', 'ASIO'}, 'Tag', 'driver', 'Units', 'normalized', 'Position', [0.7, 0.1, 0.2, 0.8],...
                'Callback', @(hObject, ~) setDriver(hObject.String{hObject.Value}));
                     
        end
        
        function playingStateSetCallback(obj, newPlayingState)
            if newPlayingState == playingStateClass('stopped')
                obj.frameSizeMenu.Visible = 'on';
                obj.deviceMenu.Visible = 'on';
                obj.driverMenu.Visible = 'on';
            else
                obj.frameSizeMenu.Visible = 'off';
                obj.deviceMenu.Visible = 'off';
                obj.driverMenu.Visible = 'off';
            end
        end
    end
    
end


