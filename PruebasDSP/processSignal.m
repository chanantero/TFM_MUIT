classdef processSignal < matlab.System & matlab.system.mixin.FiniteSource
    
    % Completar para que la salida tenga varios canales
    
    % Public, tunable properties
    properties
        
    end
    
    properties(DiscreteState)
        count
        numStoredSamples
    end
    
    properties(Nontunable)
        Fs
        frameSize % only matters if 'variable' is false
        numStoredFrames % % only matters if 'variable' is false
    end
    
    properties(Nontunable, Logical)
        variable = false
    end
    
    % Pre-computed constants
    properties(Access = private)
        storedSamples
    end
    
    methods
        function obj = processSignal(~)
            obj.Fs = 44100;
            obj.frameSize = 1024;
            obj.numStoredFrames = 1;
        end
    end
    
    methods(Access = protected)
        
        function r = stepImpl(obj, s, delay)
            
            % Calculate output r
            numChannels = size(delay, 2);
            delaySamples = round(delay*obj.Fs);
            
%             if obj.variable
                numSamples = numel(s);
%             else
%                 numSamples = obj.frameSize;
%             end
            
            indices = obj.numStoredSamples + repmat((1:numSamples)', 1, numChannels) - delaySamples;
            valid = indices > 0;
            validInd = indices(valid);
            
            availableSignal = [obj.storedSamples; s];
            
            r = zeros(numSamples, 1);
            r(valid) = availableSignal(validInd);
            
            % Liberate non-useful samples
            if obj.variable
                obj.storedSamples = availableSignal(validInd(1):end);
                obj.numStoredSamples = numel(obj.storedSamples);
            else
                obj.storedSamples = availableSignal(numSamples + 1:end); % Remove one frame and add the new one
            end
            % For when the frame size is fixed but the number of stored
            % frames is variable
            %             firstNeededToken = ceil(validInd(1)/numSamples); % Assume that validInd(j)>=validInd(i) if j > i, hence, validInd(1) == min(validInd)
            %             obj.storedSamples = availableSignal((firstNeededToken - 1)*numSamples + 1:end);
            
            % Increment counter
            obj.count = obj.count + 1;
        end
        
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            if ~obj.variable
                obj.numStoredSamples = obj.frameSize * obj.numStoredFrames;
            end
        end
        
        function resetImpl(obj)
            obj.count = 0;
            if obj.variable
                obj.numStoredSamples = 0;
                obj.storedSamples = [];
            else
                obj.storedSamples = zeros(obj.numStoredSamples, 1);
            end
        end
        
    end
end
