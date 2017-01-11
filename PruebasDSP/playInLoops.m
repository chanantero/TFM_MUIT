classdef process < matlab.System
    
    % Public, tunable properties
    properties
        increment = 1;
    end

    properties(DiscreteState)
        counter
    end

    % Pre-computed constants
    properties(Access = private)

    end

    methods 
        function obj = playInLoops(increment)
            obj.increment = increment;
        end
    end
    
    methods(Access = protected)

        function y = stepImpl(obj,u)
            y = u + obj.increment;
            obj.counter = obj.counter + 1;
        end

        function resetImpl(obj)
            obj.counter = 0;
            obj.increment = 5;
        end
    end
end
