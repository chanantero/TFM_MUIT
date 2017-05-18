classdef delayAndAttenProvider < matlab.System
    
    properties(Nontunable)
        t_change
        SampleRate
        delays
        attenuations
    end
    
    properties(DiscreteState)
        countSamples
    end
    
    properties(SetAccess = private)
        t_change_Samp
        delays_Samp
        attenuations_Samp
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
            obj.time2Samples();
        end
        
        function [delay, attenuation] = stepImpl(obj, numSamples)
            [delay, attenuation] = getDelayAndAtten(obj.countSamples, numSamples);         
            obj.countSamples = obj.countSamples + numSamples;
        end
        
        function resetImpl(obj)
            obj.countSamples = 0;
        end
        
        function releaseImpl(~)
        end
    end
    
    methods
        function obj = delayAndAttenProvider()
        end
        
        function getDelayAndAtten(obj, countSamples, numSamples)
            % countSamples: number of samples already provided
            % numSamples: number of samples of the new frame
            
            t_Samp = obj.t_change_Samp;
            del = obj.delays_Samp;
            atten = obj.attenuations_Samp;
            
            % Find out which samples of change we need
            t_ind = find(t_Samp >= (countSamples + 1) & t_Samp <= (countSamples + numSamples));
            keySamples = [t_Samp(t_ind) - countSamples; numSamples + 1];
            if keySamples(1) > 1
                keySamples = [1; keySamples];
                prevInd = find(t_Samp < countSamples + 1, 1, 'last');
                t_ind = [prevInd; t_ind];
            end
            
            % Give the delays and attenuations the convenient format for the processor
            % object
            numChann = obj.numChannels;
            delay = zeros(numSamples, numChann);
            attenuation = zeros(numSamples, numChann);
            for k = 1:numel(keySamples) - 1
                currInd = keySamples(k):keySamples(k+1)-1;
                delay(currInd, :) = repmat(del(t_ind(k), :), numel(currInd), 1);
                attenuation(currInd, :) = repmat(atten(t_ind(k), :), numel(currInd), 1);
            end
            
        end
        
        function getDelay_fullyComputed(obj, countSamples, numSamples)
            % countSamples: number of samples already provided
            % numSamples: number of samples of the new frame
            
            % Calculate samples where delay and attenuation changes
            t_Samp = ceil(obj.t_change * obj.SampleRate) + 1; % Samples were the position should change
            % Eliminate redundant information
            [t_Samp, ind] = unique(t_Samp);
            delay_Samp = obj.delays(ind, :);
            attenuation_Samp = obj.attenuations(ind, :);
            
            % Find out which samples of change we need
            t_ind = find(t_Samp >= (countSamples + 1) & t_Samp <= (countSamples + numSamp));
            keySamples = [t_Samp(t_ind) - countSamples; numSamp + 1];
            if keySamples(1) > 1
                keySamples = [1; keySamples];
                prevInd = find(t_Samp < countSamples + 1, 1, 'last');
                t_ind = [prevInd; t_ind];
            end
            
            % Give the delays and attenuations the convenient format for the processor
            % object
            numChann = obj.numChannels;
            delay = zeros(numSamples, numChann);
            attenuation = zeros(numSamples, numChann);
            for k = 1:numel(keySamples) - 1
                currInd = keySamples(k):keySamples(k+1)-1;
                delay(currInd, :) = repmat(delay_Samp(t_ind(k), :), numel(currInd), 1);
                attenuation(currInd, :) = repmat(attenuation_Samp(t_ind(k), :), numel(currInd), 1);
            end
            
        end
   
    end

    methods(Access = private)
        function time2Samples(obj)
            % Calculate samples where delay and attenuation changes
            t_Samp = ceil(obj.t_change * obj.SampleRate) + 1; % Samples were the position should change
            % Eliminate redundant information
            [obj.t_change_Samp, ind] = unique(t_Samp);
            obj.delays_Samp = obj.delays(ind, :);
            obj.attenuations_Samp = obj.attenuations(ind, :);
        end
            
    end
    
end

