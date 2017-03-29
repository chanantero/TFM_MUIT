classdef processSignal < matlab.System & matlab.system.mixin.FiniteSource
    
    % Completar para que la salida tenga varios canales
    % Completar el forward to backward delay
    
    properties(DiscreteState)
        countFrames % Number of times that the step method has been called
        countSamples % Number of samples that have been returned by the step method
        numStoredSamples % Number of samples that are currently stored
    end
    
    properties(Nontunable)
        Fs % Sampling frequency of the input signal
        delayType % Forward (real) or backward (approximation)
        numChannels % Number of channels
    end
    
    % Pre-computed constants
    properties(Access = private)
        storedSamples
        storedDelayBackward
        storedDelayForward
        storedAttenBackward
        storedAttenForward
    end
    
    methods
        function obj = processSignal(varargin)
            p = inputParser;
            
            addParameter(p, 'Fs', 44100);
            addParameter(p, 'delayType', 'backward')
            addParameter(p, 'numChannels', 1)
            
            parse(p, varargin{:})
            
            obj.Fs = p.Results.Fs;
            obj.delayType = delayTypes(p.Results.delayType);
            obj.numChannels = p.Results.numChannels;

        end
    end
    
    methods(Access = protected)
        
        function r = stepImpl(obj, s, delay, attenuation)
            % The size of the inputs delay and attenuation must have obj.numChannels columns and
            % as many rows as input samples.
            % s is a vector. The number of elements is the number of input
            % samples
            
            % Calculate output r
            numChann = obj.numChannels;
            numSamples = numel(s);
            
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
                
                storedForwAtten = obj.storedAttenForward;
                forwAtten = [storedForwAtten; attenuation];
                obj.storedAttenForward = forwAtten(end, :);
                atten = zeros(numSamples, numChann);
                
                for k = 1:numChann
                    % Delay
                    backDelay = BackFromForwDelay(storedForwDelay(:, k), delay(:, k), obj.storedDelayBackward{k}, obj.Fs);
                    delaySamples(:, k) = round(backDelay(1:numSamples)*obj.Fs);
                    obj.storedDelayBackward{k} = backDelay(numSamples+1:end);
                    
                    % Attenuation
                    backAtten = BackFromForwAtten(storedForwDelay(:, k), delay(:, k), storedForwAtten(:, k), attenuation(:,k), obj.storedAttenBackward{k}, obj.Fs);
                    atten(:, k) = backAtten(1:numSamples);
                    obj.storedAttenBackward{k} = backAtten(numSamples+1:end);
                end
                
            end
                        
            indices = obj.numStoredSamples + repmat((1:numSamples)', 1, numChann) - delaySamples;
            valid = indices > 0;
            validInd = indices(valid);
            
            r = zeros(numSamples, numChann);
            availableSignal = [obj.storedSamples; s];
            r(valid) = availableSignal(validInd);
            r = r.*atten;
            
            % Liberate non-useful samples
            obj.storedSamples = availableSignal(max(min(indices(:)), 1):end, :);
            obj.numStoredSamples = size(obj.storedSamples, 1);
            obj.countSamples = obj.countSamples + numSamples;
            
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
                obj.storedAttenForward = double.empty(0, obj.numChannels);
                obj.storedAttenBackward = cell(obj.numChannels, 1);
%             end
        end

        function releaseImpl(obj)
            % Release resources, such as file handles
            reset(obj);
        end
        
    end
end
