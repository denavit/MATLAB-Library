classdef reinf_bar < reinf
    
    properties
        x
        y
        nb = 1;
    end
    
    methods
        function obj = reinf_bar(x,y,Ab)
            obj.x = x;
            obj.y = y;
            obj.Ab = Ab;
        end
        function [x,y] = coordinates(obj)
            x = obj.x;
            y = obj.y;
        end
    end
    
end