classdef processSignal < matlab.System & matlab.system.mixin.FiniteSource
    
    properties(DiscreteState)
        countFrames % Number of times that the step method has been called
        countSamples % Number of samples that have been returned by the step method
        numStoredSamples % Number of samples that are currently stored
    end
    
    properties(Nontunable)
        Fs % Sampling frequency of the input signal
        numChannels % Number of channels
    end
    
    % Pre-computed constants
    properties(Access = private)
        storedSamples    
    end
    
    methods
        function obj = processSignal(varargin)
            p = inputParser;
            
            addParameter(p, 'Fs', 44100);
            addParameter(p, 'numChannels', 1)
            
            parse(p, varargin{:})
            
            obj.Fs = p.Results.Fs;
            obj.numChannels = p.Results.numChannels;

        end
    end
    
    methods(Access = protected)
        
        function r = stepImpl(obj, s, delays, attenuations)
            % The size of the inputs delay and attenuation must have obj.numChannels columns and
            % as many rows as input samples.
            % s is a vector. The number of elements is the number of input
            % samples
            
            % Calculate output r
            numChann = obj.numChannels;
            numSamples = numel(s);
            
            availableSignal = [obj.storedSamples; s];
            
            numComp = size(delays, 3);
            
            r = zeros(numSamples, numChann, numComp);
            minIndices = zeros(numComp, 1);
            for l = 1:numComp
                
                delay = delays(:, :, l);
                delaySamples = round(delay*obj.Fs);
                atten = attenuations(:, :, l);
                               
                % Apply Delay
                indices = obj.numStoredSamples + repmat((1:numSamples)', 1, numChann) - delaySamples;
                minIndices(l) = min(indices(:));
                valid = indices > 0;
                
                r_comp = zeros(numSamples, numChann);
                r_comp(valid) = availableSignal(indices(valid));
                
                % Apply Attenuation
                r_comp = r_comp.*atten;
                
                r(:,:,l) = r_comp;
                
            end
            r = sum(r, 3);
            
            
            % Liberate non-useful samples
            obj.storedSamples = availableSignal(max(min(minIndices), 1):end, :);
            
            % Update discrete state
            obj.numStoredSamples = size(obj.storedSamples, 1);
            obj.countSamples = obj.countSamples + numSamples;
            obj.countFrames = obj.countFrames + 1;
        end
        
        function resetImpl(obj)
            obj.countFrames = 0;
            obj.countSamples = 0;
            obj.numStoredSamples = 0;
            obj.storedSamples = [];           
        end

        function releaseImpl(obj)
            % Release resources, such as file handles
            reset(obj);
        end
        
    end
end
