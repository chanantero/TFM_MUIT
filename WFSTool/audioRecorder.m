classdef audioRecorder < matlab.System
    % Untitled Add summary here
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.

    % Public, tunable properties
    properties
        Driver
        Device
        SamplesPerFrame
        SampleRate
        channelMapping
    end

    properties(DiscreteState)
        count
    end
    
    properties(Dependent)
        NumChannels
    end

    properties(SetAccess = private)
        audioDeviceReaderObj
        recorded
    end
    
    % Getters and setters
    methods
        function NumChannels = get.NumChannels(obj)
            NumChannels = obj.audioDeviceReaderObj.NumChannels;
        end
    end

    methods(Access = protected)
        function setupImpl(obj)
           obj.audioDeviceReaderObj.Driver = obj.Driver;
           obj.audioDeviceReaderObj.Device = obj.Device;
           obj.audioDeviceReaderObj.SamplesPerFrame = obj.SamplesPerFrame;
           obj.audioDeviceReaderObj.SampleRate = obj.SampleRate;
           obj.audioDeviceReaderObj.ChannelMappingSource = 'Property';
           obj.audioDeviceReaderObj.ChannelMapping = obj.channelMapping;
           
           
           setup(obj.audioDeviceReaderObj);
        end

        function stepImpl(obj)
            obj.recorded = [obj.recorded; step(obj.audioDeviceReaderObj)];
            
            obj.count = obj.count + 1;
        end

        function resetImpl(obj)
            reset(obj.audioDeviceReaderObj);
            obj.recorded = zeros(0, obj.NumChannels);
            obj.count = 0;
        end
        
        function releaseImpl(obj)
            release(obj.audioDeviceReaderObj);
        end
    end
    
    methods
        function obj = audioRecorder()
            obj.audioDeviceReaderObj = audioDeviceReader;
            
            obj.setDefaultProperties();
        end
        
        function devices = getAudioDevices(obj)
%             devices = getAudioDevices(obj.audioDeviceReaderObj);
            
            % If there is not set method for driver, use the next two lines
            aux = audioDeviceReader('Driver', obj.Driver);
            devices = getAudioDevices(aux);
        end
    end
    
    methods(Access = private)
        function setDefaultProperties(obj)
            obj.Driver = 'DirectSound';
            obj.Device = 'Default';
            obj.SampleRate = 44100;
            obj.SamplesPerFrame = 44100;
        end
        
    end
end
