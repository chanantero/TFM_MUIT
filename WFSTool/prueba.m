classdef prueba < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        var2
    end
    
    properties(Dependent)
        var1
    end
    
    methods
        function obj = prueba()
            obj.var1 = 1:4;
        end
        
        function a = get.var1(obj)
            a = cell2mat(obj.var2) + 1;
        end
        
        function set.var1(obj, value)
            obj.var2 = num2cell(value);
        end
        
    end
    
end

