classdef delayAndAttenProvider < matlab.System
    
    properties(Nontunable)
        mode
        FrameSize
        % mode == realTime
        getDelayFun % Cell array of functions with as many elements as processors
        getAttenFun % Cell array of functions with as many elements as processors
        % mode == predefined
        SampleRate
        t_change
        delays
        attenuations
    end
    
    properties(DiscreteState)
        count
    end
    
    properties(SetAccess = private)
        t_change_Samp
        delays_Samp
        attenuations_Samp
        getDelayAndAttenFun
    end
    
    
    properties(Dependent)
        numChannels
    end
    
    % Getters and setters
    methods
        function numChannels = get.numChannels(obj)
            numChannels = size(obj.delays, 2);
        end
    end
    
    methods(Access = protected)
        function setupImpl(obj)
            obj.time2Samples_update();
            if obj.mode == timeInteractionTypes('realTime');
                obj.getDelayFun() = @() obj.stepRealTime();
            else
                obj.getDelayFun() = @() obj.stepPredefined();
            end
        end
        
        function [delay, attenuation] = stepImpl(obj)
            [delay, attenuation] = obj.getDelayAndAttenFun();
            obj.count = obj.count + 1;
        end
        
        function resetImpl(obj)
            obj.countSamples = 0;
        end
        
        function releaseImpl(~)
        end
    end
    
    methods
        function obj = delayAndAttenProvider()
            obj.mode = timeInteractionTypes('realTime');
        end
        
        function [delay, attenuation] = stepRealTime(obj)
            delay = repmat(permute(obj.getDelayFun(), [3, 1, 2]), obj.FrameSize, 1);
            attenuation = repmat(permute(obj.getAttenFun(), [3, 1, 2]), obj.FrameSize, 1);
        end
        
        function [delay, attenuation] = stepPredefined(obj)
            firstSample = obj.count * obj.FrameSize + 1;
            lastSample = firstSample + obj.FrameSize - 1;
            
            [delay, attenuation] = getDelayAndAtten(firstSample, lastSample);
        end
        
        function [delay, attenuation] = getDelayAndAtten(obj, firstSample, lastSample)
            % countSamples: number of samples already provided
            % numSamples: number of samples of the new frame
            
            t_Samp = obj.t_change_Samp;
            del = obj.delays_Samp;
            atten = obj.attenuations_Samp;
            
            % Find out which samples of change we need
            [keySamples, t_ind] = vectorSubSelection(t_Samp, firstSample, lastSample);
            
            % If it is empty or the first index is bigger than 1, we need
            % to add a first keySample of change
            if isempty(keySamples) || keySamples(1) > 1
                keySamples = [1; keySamples];
                prevInd = find(t_Samp < firstSample, 1, 'last');
                % If prevInd is empty, it is because there are not previous
                % change samples. We will assume the value it is the same
                % as the first
                if isempty(prevInd)
                    prevInd = t_ind(1);
                end
                t_ind = [prevInd; t_ind];
            end
            
            % Give the delays and attenuations the convenient format for the processor
            % object
            valueInd = changeIndices2Vector([keySamples; lastSample - firstSample + 1], t_ind);
            delay = del(valueInd, :);
            attenuation = atten(valueInd, :);
            
        end
        
        function [delay, attenuation] = getDelay_fullyComputed(obj, firstSample, lastSample)
            % countSamples: number of samples already provided
            % numSamples: number of samples of the new frame
            
            % Calculate samples where delay and attenuation changes
            [t_Samp, delay_Samp, attenuation_Samp] = delayAndAttenProvider.time2Samples(obj.SampleRate, obj.t_change, obj.delays, obj.attenuations);
            
            % Find out which samples of change we need
            [keySamples, t_ind] = vectorSubSelection(t_Samp, firstSample, lastSample);
            
            % If it is empty or the first index is bigger than 1, we need
            % to add a first keySample of change
            if isempty(keySamples) || keySamples(1) > 1
                keySamples = [1; keySamples];
                prevInd = find(t_Samp < firstSample, 1, 'last');
                % If prevInd is empty, it is because there are not previous
                % change samples. We will assume the value it is the same
                % as the first
                if isempty(prevInd)
                    prevInd = t_ind(1);
                end
                t_ind = [prevInd; t_ind];
            end
            
            % Give the delays and attenuations the convenient format for the processor
            % object
            valueInd = changeIndices2Vector([keySamples; lastSample - firstSample + 1], t_ind);
            delay = delay_Samp(valueInd, :);
            attenuation = attenuation_Samp(valueInd, :);
            
        end
        
    end
    
    methods(Access = private)
        function time2Samples_update(obj)
            [obj.t_Samp, obj.delays_Samp, obj.attenuations_Samp] = delayAndAttenProvider.time2Samples(obj.SampleRate, obj.t_change, obj.delays, obj.attenuations);
        end
        
    end
    
    
    methods(Static)
        function [timeSamp, varargout] = time2Samples(SampleRate, time, varargin)
            
            % Calculate samples where delay and attenuation changes
            timeSamp = ceil(time * SampleRate) + 1; % Samples were the position should change
            % Eliminate redundant information
            [timeSamp, ind] = unique(timeSamp);
            
            numVal = nargin - 2;
            varargout = cell(numVal, 1);
            for k = 1:numVal
                varargout{k} = varargin{k}(ind, :);
            end
            
        end
        
        function [selInd, valueInd] = vectorSubSelection(x, firstSample, lastSample)
            % x must be an sorted vector in ascending order
            valueInd = find(x >= firstSample & t_Samp <= lastSample);
            selInd = x(valueInd) - firstSample + 1;
        end
        
        function values = changeIndices2Vector(changeIndices, valueIndices)
            % changeIndices must have one more element than valueIndices.
            % That last element indicates the ending
            
            for k = 1:numel(changeIndices) - 1
                currInd = changeIndices(k):changeIndices(k+1) - 1;
                values = repmat(valueIndices(k, :), numel(currInd), 1);
            end
            
        end
    end
end


