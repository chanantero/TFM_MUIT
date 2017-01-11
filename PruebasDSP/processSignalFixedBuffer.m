classdef processSignalFixedBuffer < matlab.System & matlab.system.mixin.FiniteSource
    
    % Public, tunable properties
    properties
        
    end

    properties(DiscreteState)
        count
    end

    properties(Nontunable)
        Fs
        frameSize
        numStoredFrames
    end
    
    % Pre-computed constants
    properties(Access = private)
        storedSamples
        numStoredSamples
    end
    
    methods
        function obj = processSignalFixedBuffer(~)
            obj.Fs = 44100;
            obj.frameSize = 1024;
            obj.numStoredFrames = 1;
        end
    end
    
    methods(Access = protected)     
        
        function r = stepImpl(obj, s, delay)
    
            % Calculate output r
            delaySamples = round(delay*obj.Fs);
            
            numSamples = obj.frameSize;
            
            indices = obj.numStoredSamples + (1:numSamples)' - delaySamples;
            valid = indices > 0;
            validInd = indices(valid);
            
            availableSignal = [obj.storedSamples; s];
            
            r = zeros(numSamples, 1);
            r(valid) = availableSignal(validInd);
            
            % Liberate non-useful samples
            % For when the frame size is fixed but the number of stored
            % frames is variable
%             firstNeededToken = ceil(validInd(1)/numSamples); % Assume that validInd(j)>=validInd(i) if j > i, hence, validInd(1) == min(validInd)
%             obj.storedSamples = availableSignal((firstNeededToken - 1)*numSamples + 1:end);
            obj.storedSamples = availableSignal(numSamples + 1:end); % Remove one frame and add the new one
            
            % Increment counter
            obj.count = obj.count + 1;
        end

        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            obj.numStoredSamples = obj.frameSize * obj.numStoredFrames;
        end

        function resetImpl(obj)
            obj.count = 0;
            obj.storedSamples = zeros(obj.numStoredSamples, 1);
        end

    end
end
