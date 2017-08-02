classdef pruebaDependent < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Dependent)
        uno
    end
    
    properties
        dos
    end
    
    methods
        function obj = pruebaDependent()
            obj.dos = 4;
        end
        
        function uno = get.uno(obj)
            uno = obj.dos;
        end
        
        function set.uno(obj, value)
            obj.dos = value;
        end
    end
    
end

