classdef audioPlayer < matlab.System
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
        function obj = audioPlayer(varargin)   
            p = inputParser;
            
            addParameter(p, 'Fs', 44100);
            addParameter(p, 'frameSize', 1024);
            addParameter(p, 'numChannels', 2);
            
            parse(p, varargin{:});
            
            obj.Fs = p.Results.Fs;
            obj.frameSize = p.Results.frameSize;
            obj.numChannels = p.Results.numChannels;
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
            frSize = obj.frameSize;
            numFrames = floor(size(stored, 1)/frSize);
            if numFrames > 0 % Send it for the reproduction
                for k = 1:numFrames
                    sampleIndices = (k-1)*frSize+1:k*frSize;
                    frame = stored(sampleIndices, :);
                    play(obj.deviceWriter, frame);
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
end
