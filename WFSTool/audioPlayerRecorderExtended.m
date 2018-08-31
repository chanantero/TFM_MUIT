classdef audioPlayerRecorderExtended < matlab.System
    % When the stored data is big enough, it sends it to be reproduced.

    % Public, tunable properties
    properties

    end

    properties(Nontunable)
        SampleRate
        SamplesPerFrame
        PlayerChannelMapping
        RecorderChannelMapping
        Device
        
        % Signal provider configuration
        mode % originType('file'), originType('custom'), originType('func'), originType('sinusoidal')
            % Variables that are exclusive for the mode == file
        filename        
            % Variables that are exclusive for the mode == custom
        customSignal
            % Variables that are exclusive for the mode == func
        signalFunc
    end
    
    properties(DiscreteState)
        count
    end

    properties(Dependent)
        numPlayerChannels
        numRecorderChannels
        numSignalChannels
        numSamples
        DefaultPlayerChannels
        DefaultRecorderChannels
    end
    
    % Pre-computed constants
    properties(SetAccess = private)
        recSignal
        playerRecorder
        signProv
    end

    methods
        function obj = audioPlayerRecorderExtended(varargin)   
            p = inputParser;
            
            addParameter(p, 'SampleRate', 44100);
            addParameter(p, 'SamplesPerFrame', 1024);
            addParameter(p, 'PlayerChannelMapping', []);
            addParameter(p, 'RecorderChannelMapping', []);
            addParameter(p, 'Device', 'Default');

            parse(p, varargin{:});
            
            obj.SampleRate = p.Results.SampleRate;
            obj.SamplesPerFrame = p.Results.SamplesPerFrame;
            obj.PlayerChannelMapping = p.Results.PlayerChannelMapping;
            obj.RecorderChannelMapping = p.Results.RecorderChannelMapping;
            obj.Device = p.Results.Device;
            obj.mode = originType('file');
            obj.filename = '';
            obj.customSignal = [];

            obj.playerRecorder = audioPlayerRecorder;
            obj.signProv = signalProvider;
        end
        
        function playAndRecord(obj)
            setup(obj);
            reset(obj);
            
            step(obj.playerRecorder, zeros(obj.SamplesPerFrame, obj.numSignalChannels));
            
            numFrames = ceil(obj.numSamples/obj.SamplesPerFrame);
            for fr = 1:numFrames
                step(obj);
            end
            
            release(obj);
        end
        
        function devices = getAudioDevices(obj)
            devices = getAudioDevices(obj.playerRecorder);
        end
        
        function numSamples = get.numSamples(obj)
            numSamples = obj.signProv.numSamples;
        end
        
        function numPlayerChannels = get.numPlayerChannels(obj)
            if obj.DefaultPlayerChannels
                inf = info(obj.playerRecorder);
                numPlayerChannels = inf.MaximumPlayerChannels;
            else
                numPlayerChannels = numel(obj.PlayerChannelMapping);
            end
        end
        
        function numRecorderChannels = get.numRecorderChannels(obj)
            if obj.DefaultRecorderChannels
                inf = info(obj.playerRecorder);
                numRecorderChannels = inf.MaximumRecorderChannels;
            else
                numRecorderChannels = numel(obj.RecorderChannelMapping);
            end
        end
        
        function numSignalChannels = get.numSignalChannels(obj)
            numSignalChannels = obj.signProv.numChannels;
        end
        
        function DefaultPlayerChannels = get.DefaultPlayerChannels(obj)
            DefaultPlayerChannels = isempty(obj.PlayerChannelMapping);
        end
        
        function DefaultRecorderChannels = get.DefaultRecorderChannels(obj)
            DefaultRecorderChannels = isempty(obj.RecorderChannelMapping);
        end
        
    end
    
    methods(Access = protected)
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            obj.setAudioPlayerRecorderProperties();
            obj.setSignalProviderProperties();
            
            obj.recSignal = zeros(obj.signProv.numSamples, obj.numRecorderChannels);
            
            % Don't do the setup now if it's not necessary. Let's avoid
            % possible sources of desynchronization
%             setup(obj.playerRecorder, zeros(obj.SamplesPerFrame, obj.numSignalChannels));
%             setup(obj.signProv);
        end

        function [numUnderrun, numOverrun] = stepImpl(obj)
            x = step(obj.signProv);
            
            startSample = obj.count * obj.SamplesPerFrame + 1;
            endSample = startSample + obj.SamplesPerFrame - 1;
            ind = startSample:endSample;
            [obj.recSignal(ind, :), numUnderrun, numOverrun] = step(obj.playerRecorder, x); 
            
            obj.count = obj.count + 1;
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            obj.count = 0;
            obj.recSignal = double.empty(0, obj.numRecorderChannels);
        end

        function releaseImpl(obj)
            % Release resources, such as file handles
            release(obj.playerRecorder);
            release(obj.signProv);
        end
                
    end
    
    methods(Access = private)
        function setAudioPlayerRecorderProperties(obj)
            obj.playerRecorder.SampleRate = obj.SampleRate;
            obj.playerRecorder.Device = obj.Device;
            
            obj.playerRecorder.PlayerChannelMapping = obj.PlayerChannelMapping;
            obj.playerRecorder.RecorderChannelMapping = obj.RecorderChannelMapping;
        end
        
        function setSignalProviderProperties(obj)
            obj.signProv.mode = obj.mode;
            obj.signProv.SamplesPerFrame = obj.SamplesPerFrame;
            obj.signProv.SampleRate = obj.SampleRate;
            obj.signProv.FileName = obj.filename;
            obj.signProv.customSignal = obj.customSignal;
        end
    end
end
