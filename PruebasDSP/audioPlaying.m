classdef audioPlaying < matlab.System
    % When the stored data is big enough, it sends it to be reproduced.

    % Public, tunable properties
    properties

    end

    properties(Nontunable)
        Fs
        frameSize
        numChannels
    end

    properties(DiscreteState)
        count
    end

    % Pre-computed constants
    properties(Access = private)
        storedSamples
        deviceWriter
    end

    methods
        function obj = audioPlaying(~)
            obj.Fs = 44100;
            obj.frameSize = 1024;
            obj.numChannels = 2;
        end
    end
    
    methods(Access = protected)
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            obj.deviceWriter = audioDeviceWriter('SampleRate', obj.Fs);
            setup(obj.deviceWriter, zeros(obj.frameSize, obj.numChannels));
        end

        function stepImpl(obj, signal)
            stored = [obj.storedSamples; signal];
            obj.storedSamples = stored;
            frSize = obj.frameSize;
            numFrames = floor(numel(stored)/frSize);
            if numFrames > 0 % Send it for the reproduction
                for k = 1:numFrames
                    sampleIndices = (k-1)*frSize+1:k*frSize;
                    frame = stored(sampleIndices, :);
                    play(obj.deviceWriter, frame);
                    obj.count = obj.count + 1;
                end
            end            
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            obj.count = 0;
        end

        function releaseImpl(obj)
            % Release resources, such as file handles
            release(deviceWriter);
        end
        
        
    end
end
