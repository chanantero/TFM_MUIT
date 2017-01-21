classdef processSignal < matlab.System & matlab.system.mixin.FiniteSource
    
    % Completar para que la salida tenga varios canales
    % Completar el forward to backward delay
    
    % Public, tunable properties
    properties
        
    end
    
    properties(DiscreteState)
        countFrames
        countSamples
        numStoredSamples
    end
    
    properties(Nontunable)
        Fs
        frameSize % only matters if 'variable' is false
        numStoredFrames % % only matters if 'variable' is false
        delayType
        numChannels
    end
    
    properties(Nontunable, Logical)
        variable = false
    end
    
    % Pre-computed constants
    properties(Access = private)
        storedSamples
        storedDelayBackward
        storedDelayForward
    end
    
    methods
        function obj = processSignal(varargin)
            p = inputParser;
            
            addParameter(p, 'Fs', 44100);
            addParameter(p, 'frameSize', 1024);
            addParameter(p, 'numStoredFrames', 1);
            addParameter(p, 'variable', false);
            addParameter(p, 'delayType', 'backward')
            addParameter(p, 'numChannels', 1)
            
            parse(p, varargin{:})
            
            obj.Fs = p.Results.Fs;
            obj.frameSize = p.Results.frameSize;
            obj.numStoredFrames = p.Results.numStoredFrames;
            obj.variable = p.Results.variable;
            obj.delayType = delayTypes(p.Results.delayType);
            obj.numChannels = p.Results.numChannels;

        end
    end
    
    methods(Access = protected)
        
        function r = stepImpl(obj, s, delay)
            
            % Calculate output r
            numChann = obj.numChannels;
            
%             if obj.variable
                numSamples = numel(s);
%             else
%                 numSamples = obj.frameSize;
%             end
            
            if obj.delayType % delay is the backward delay
                delaySamples = round(delay*obj.Fs);                
            else % delay is the forward delay
                % Calculate the backward delay
                % Each time a frame arrives, the backward delay is
                % calculated and appended to the already stored backward
                % delay.
                % storedDelayBackward must be initialized to cell array
                % storedDelayForward must be initialized to empty array
                storedForwDelay = obj.storedDelayForward;
                forwDelay = [storedForwDelay; delay];
                obj.storedDelayForward = forwDelay(end, :);
                delaySamples = zeros(numSamples, numChann);
                
                for k = 1:numChann
                    backDelay = BackFromForwDelay(storedForwDelay(:, k), delay(:, k), obj.storedDelayBackward{k}, obj.Fs);
                    delaySamples(:, k) = round(backDelay(1:numSamples)*obj.Fs);
                    obj.storedDelayBackward{k} = backDelay(numSamples+1:end);
                end
                
            end
            
            indices = obj.numStoredSamples + repmat((1:numSamples)', 1, numChann) - delaySamples;
            valid = indices > 0;
            validInd = indices(valid);
            
            availableSignal = [obj.storedSamples; s];
            
            r = zeros(numSamples, numChann);
            r(valid) = availableSignal(validInd);
            
            % Liberate non-useful samples
            if obj.variable
                obj.storedSamples = availableSignal(max(min(indices), 1):end, :);
                obj.numStoredSamples = size(obj.storedSamples, 1);
                obj.countSamples = obj.countSamples + numSamples;
            else
                obj.storedSamples = availableSignal(numSamples + 1:end, :); % Remove one frame and add the new one
                obj.countSamples = obj.countSamples + obj.frameSize;
            end
            % For when the frame size is fixed but the number of stored
            % frames is variable
            %             firstNeededToken = ceil(validInd(1)/numSamples); % Assume that validInd(j)>=validInd(i) if j > i, hence, validInd(1) == min(validInd)
            %             obj.storedSamples = availableSignal((firstNeededToken - 1)*numSamples + 1:end);
            
            % Increment counter
            obj.countFrames = obj.countFrames + 1;
        end
        
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            if ~obj.variable
                obj.numStoredSamples = obj.frameSize * obj.numStoredFrames;
            end
        end
        
        function resetImpl(obj)
            obj.countFrames = 0;
            obj.countSamples = 0;
            if obj.variable
                obj.numStoredSamples = 0;
                obj.storedSamples = [];
            else
                obj.storedSamples = zeros(obj.numStoredSamples, obj.numChannels);
            end
            
%             if ~obj.delayType
                obj.storedDelayForward = double.empty(0, obj.numChannels);
                obj.storedDelayBackward = cell(obj.numChannels, 1);
%             end
        end

        function releaseImpl(obj)
            % Release resources, such as file handles
            reset(obj);
        end
        
    end
end
