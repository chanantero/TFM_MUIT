classdef audioPlayer < matlab.System
    % When the stored data is big enough, it sends it to be reproduced.

    properties(Nontunable, Logical)
        DefaultChannels = true;
    end

    % Public, tunable properties
    properties

    end

    properties(Nontunable)
        Fs
        frameSize
        channelMapping
        device
        driver        
    end
    
    properties(DiscreteState)
        count
    end

    properties(Dependent)
        numChannels
    end
    
    % Pre-computed constants
    properties(SetAccess = private)
        storedSamples
        deviceWriter
    end

    methods
        function obj = audioPlayer(varargin)   
            p = inputParser;
            
            addParameter(p, 'Fs', 44100);
            addParameter(p, 'frameSize', 1024);
            addParameter(p, 'channelMapping', 2);
            addParameter(p, 'device', 'Default');
            addParameter(p, 'driver', 'DirectSound');
            addParameter(p, 'DefaultChannels', true);

            parse(p, varargin{:});
            
            obj.Fs = p.Results.Fs;
            obj.frameSize = p.Results.frameSize;
            obj.channelMapping = p.Results.channelMapping;
            obj.device = p.Results.device;
            obj.driver = p.Results.driver;
            obj.DefaultChannels = p.Results.DefaultChannels;

            obj.deviceWriter = audioDeviceWriter;
        end
        
%         function devices = getAvailableDevices(obj)
% %             devices = getAudioDevices(obj.deviceWriter);
%             
%             % If there is not set method for driver, use the next two lines
%             aux = audioDeviceWriter('Driver', obj.driver);
%             devices = getAudioDevices(aux);
%         end
        
        function devices = getAudioDevices(obj)
            devices = getAudioDevices(obj.deviceWriter);
            
%             % If there is not set method for driver, use the next two lines
%             aux = audioDeviceWriter('Driver', obj.driver);
%             devices = getAudioDevices(aux);
        end
        
        function numChannels = get.numChannels(obj)
            if obj.DefaultChannels
                inf = info(obj.deviceWriter);
                numChannels = inf.MaximumOutputChannels;
            else
                numChannels = numel(obj.channelMapping);
            end
        end
        
    end
    
    methods(Access = protected)
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            obj.setDeviceWriterProperties();
            
            setup(obj.deviceWriter, zeros(obj.frameSize, obj.numChannels));
        end

        function [numUnderrun, t] = stepImpl(obj, signal)
            stored = [obj.storedSamples; signal];
            frSize = obj.frameSize;
            numFrames = floor(size(stored, 1)/frSize);
            inf = info(obj.deviceWriter);
            maxNumOutputChann = inf.MaximumOutputChannels;
            if numFrames > 0 % Send it for the reproduction
                for k = 1:numFrames
                    sampleIndices = (k-1)*frSize+1:k*frSize;
                    frame = stored(sampleIndices, 1:maxNumOutputChann);
                    numUnderrun = play(obj.deviceWriter, frame);
                    t = tic;
                    obj.count = obj.count + 1;
                end
                stored = stored(sampleIndices(end)+1:end, :);
            end            
            obj.storedSamples = stored;
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            obj.count = 0;
            obj.storedSamples = double.empty(0, obj.numChannels);
        end

        function releaseImpl(obj)
            % Release resources, such as file handles
            release(obj.deviceWriter);
        end
                
    end
    
    methods(Access = private)
        function setDeviceWriterProperties(obj)
            obj.deviceWriter.SampleRate = obj.Fs;
            obj.deviceWriter.Driver = obj.driver;
            obj.deviceWriter.Device = obj.device;
            
            if obj.DefaultChannels
                obj.deviceWriter.ChannelMappingSource = 'Auto';
            else
                obj.deviceWriter.ChannelMappingSource = 'Property';
                obj.deviceWriter.ChannelMapping = obj.channelMapping;
            end

            
            %             obj.deviceWriter.SupportVariableSizeInput = true; % Por qué
%             no iba sin esto y ahora sí que va?
        end
    end
end
