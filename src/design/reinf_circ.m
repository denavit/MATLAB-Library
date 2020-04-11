classdef reinf_circ < reinf
    
    properties
        xc      % x coordiante of the center of bar pattern
        yc      % y coordiante of the center of bar pattern
        rc      % radius from center of the bar patter to the center of the bar
        nb      % number of bars around the circumference 
    end
    
    methods
        function obj = reinf_circ(xc,yc,rc,nb,Ab,db)
            obj.xc = xc;
            obj.yc = yc;
            obj.rc = rc;
            obj.nb = nb;
            obj.Ab = Ab;
            if nargin > 5
                obj.db = db;
            end
        end
        function [x,y] = coordinates(obj)
            angles = linspace(0,2*pi,obj.nb+1);
            angles = angles(1:(end-1));
            x = obj.xc + cos(angles)*obj.rc;
            y = obj.yc + sin(angles)*obj.rc;
        end
        function n = num_bars(obj)
            n = obj.nb;
        end
    end
end