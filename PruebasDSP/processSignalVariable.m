classdef processSignalVariable < matlab.System & matlab.system.mixin.FiniteSource
    
    % Public, tunable properties
    properties
        
    end

    properties(DiscreteState)
        count
        numStoredSamples
    end

    properties(Nontunable)
        Fs
    end
    
    % Pre-computed constants
    properties(Access = private)
        storedSamples
    end
    
    methods
        function obj = processSignalVariable(~)
            obj.Fs = 44100;
        end
    end
    
    methods(Access = protected)
        
        function r = stepImpl(obj, s, delay)
    
            % Calculate output r
            delaySamples = round(delay*obj.Fs);
            
            numSamples = numel(s);
            indices = obj.numStoredSamples + (1:numSamples)' - delaySamples;
            valid = indices > 0;
            validInd = indices(valid);
            
            availableSignal = [obj.storedSamples; s];
            
            r = zeros(numSamples, 1);
            r(valid) = availableSignal(validInd);
            
            % Liberate non-useful samples
            obj.storedSamples = availableSignal(min(validInd):end);
            obj.numStoredSamples = numel(obj.storedSamples);
            
            % Increment counter
            obj.count = obj.count + 1;
        end

        function resetImpl(obj)
            obj.count = 0;
            obj.numStoredSamples = 0;
            obj.storedSamples = [];
        end

    end
end
